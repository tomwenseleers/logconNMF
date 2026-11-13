sinmf_core_count <- function(ncores) {
  if(is.null(ncores)){ncores<-parallel::detectCores(logical=FALSE);if(is.na(ncores))ncores<-1L}
  if(!is.numeric(ncores)||length(ncores)!=1L||!is.finite(ncores)||ncores<1||
      ncores!=floor(ncores)||ncores>.Machine$integer.max)stop("ncores must be a positive integer.")
  as.integer(ncores)
}

sinmf_runtime <- function(ncores,n,workers=NULL) {
  ncores<-sinmf_core_count(ncores)
  if(is.null(workers))workers<-min(n,ncores)
  workers<-min(n,ncores,sinmf_core_count(workers));threads<-ncores%/%workers
  env<-new.env(parent=emptyenv());env$cl<-NULL
  env$budget<-list(ncores=ncores,workers=workers,threads=threads)
  env$blas<-RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(1L)
  success<-FALSE
  on.exit(if(!success)RhpcBLASctl::blas_set_num_threads(env$blas),add=TRUE)
  if(workers>1L) {
    keys<-c("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","OMP_WAIT_POLICY","R_LIBS")
    old<-Sys.getenv(keys,unset=NA_character_)
    on.exit({Sys.unsetenv(keys[is.na(old)]);do.call(Sys.setenv,as.list(old[!is.na(old)]))},add=TRUE)
    Sys.setenv(OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,OMP_NUM_THREADS=threads,OMP_WAIT_POLICY="PASSIVE",
      R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
    env$cl<-parallel::makePSOCKcluster(workers,rscript_args="--vanilla",useXDR=FALSE)
    ok<-FALSE;on.exit(if(!ok)parallel::stopCluster(env$cl),add=TRUE)
    initialize<-function(paths){.libPaths(paths);loadNamespace("logconNMF");NULL}
    environment(initialize)<-baseenv()
    parallel::clusterCall(env$cl,initialize,.libPaths())
    ok<-TRUE
  }
  success<-TRUE;env
}
sinmf_close_runtime <- function(runtime) {
  on.exit(RhpcBLASctl::blas_set_num_threads(runtime$blas),add=TRUE)
  if(!is.null(runtime$cl))parallel::stopCluster(runtime$cl)
}
sinmf_apply <- function(runtime,tasks,fun) {
  if(is.null(runtime)||is.null(runtime$cl))lapply(tasks,fun)else
    parallel::parLapplyLB(runtime$cl,tasks,fun,chunk.size=1L)
}

# A task contains one sample's compact spectrum slice and paths. It never
# captures the state, shared observation cache, or every sample's design.
sinmf_tasks <- function(state,samples,operation,roots=NULL,...) {
  lapply(samples,function(s) {
    f<-readRDS(state$paths[s]);ids<-state$map[f$gid]
    list(path=state$paths[s],directory=state$store$manifest$directory[s],
      P=state$P[ids,,drop=FALSE],scale=state$scale[ids],operation=operation,
      selected=if(is.null(roots))seq_along(ids)else which(ids%in%roots),
      threads=if(is.null(state$runtime))f$ncores else state$runtime$budget$threads,
      groups=state$channel_groups,
      aggregate=if(operation!="profile"||is.null(state$channel_groups)||!f$kind%in%1:2||length(f$phi)!=1L||
        !is.null(dim(f$weights)))NULL else sinmf_sample_data(state$store,s,aggregate=TRUE),
      extra=list(...))
  })
}

sinmf_sample_operation <- function(task) {
  f<-readRDS(task$path);old_threads<-f$ncores;gid<-f$gid
  f$ncores<-as.integer(task$threads);f$gid<-seq_along(gid)
  meta<-readRDS(file.path(task$directory,"metadata.rds"))
  tiny<-list(map=seq_along(gid),P=task$P,scale=task$scale)
  if(task$operation%in%c("score","gradient","background"))
    return(sinmf_stream_operation(task,f,meta,tiny,old_threads,gid))
  if(task$operation=="profile"&&!is.null(task$aggregate))
    return(sinmf_aggregate_profile(task,f,meta,tiny,old_threads,gid))
  Y<-sinmf_read_matrix_file(file.path(task$directory,"observed.bin"),meta$n,meta$m)
  mu<-sinmf_mean(f,tiny)
  if(task$operation=="profile") {
    block<-task$selected
    for(k in block[order(f$area[block]*task$scale[block],decreasing=TRUE)]) {
      at<-which(f$H[,k]>0);if(length(at)<3L)next
      h<-as.numeric(f$H[at,k]);sp<-task$P[k,]*f$area[k]*task$scale[k]
      W<-sinmf_observation_block(f$weights,at)
      z<-slcn_log_profile(Y[at,,drop=FALSE],mu[at,,drop=FALSE],h,sp,
        f$kind,f$phi,f$theta,as.numeric(W),length(at),censor=f$censor,
        maxit=as.integer(task$extra$steps),q=as.integer(task$extra$q),nthreads=f$ncores,mode=2L,
        time=f$aligned_time[at],exposure=f$exposure[at])
      hnew<-z$profile;scale<-sum(hnew*f$quadrature[at]);if(!is.finite(scale)||scale<=0)next
      nextmu<-mu[at,,drop=FALSE]+(f$exposure[at]*(hnew-h))%o%sp
      delta<-slcn_loglik(Y[at,,drop=FALSE],nextmu,f$kind,f$phi,f$theta,as.numeric(W),f$ncores,f$censor)-
        slcn_loglik(Y[at,,drop=FALSE],mu[at,,drop=FALSE],f$kind,f$phi,f$theta,as.numeric(W),f$ncores,f$censor)
      if(is.finite(delta)&&delta>= -1e-7){f$H[at,k]<-hnew/scale;f$area[k]<-f$area[k]*scale;mu[at,]<-nextmu}
    }
  } else if(task$operation=="noise"&&f$kind==3L&&isTRUE(f$spec$theta_update)) {
    for(iteration in 1:3) {
      score<-slcn_nb_score(Y,mu,f$theta,as.numeric(f$weights),f$ncores)
      move<-if(is.finite(score[2L])&&score[2L]<0)-score[1L]/score[2L]else sign(score[1L])
      move<-max(-1,min(1,move));theta0<-f$theta;reference<-sinmf_loglik(f,Y,mu)
      for(half in 0:12) {
        f$theta<-max(1e-6,min(1e8,theta0*exp(move/2^half)))
        if(sinmf_loglik(f,Y,mu)>=reference)break
        f$theta<-theta0
      }
      if(abs(log(f$theta/theta0))<1e-4)break
    }
    f$spec$theta<-f$theta
  }
  f$gid<-gid;f$ncores<-old_threads
  sinmf_checkpoint_sample(f,task$directory,task$operation)
}

sinmf_parallel_pass <- function(state,samples,operation,roots=NULL,...) {
  result<-sinmf_apply(state$runtime,sinmf_tasks(state,samples,operation,roots,...),sinmf_sample_operation)
  state$paths[samples]<-unlist(result,use.names=FALSE);state
}


sinmf_slice_observation <- function(x,channels)if(is.null(dim(x)))x else x[,channels,drop=FALSE]
sinmf_slice_noise <- function(x,channels)if(length(x)==1L)x else x[channels]
sinmf_block_mean <- function(f,P,scale,channels,rows=seq_len(nrow(f$H))) {
  S<-P[,channels,drop=FALSE]*(f$area*scale)
  as.matrix(f$H[rows,,drop=FALSE]%*%S)*f$exposure[rows]+
    f$background[rows,,drop=FALSE]%*%f$alpha[,channels,drop=FALSE]*f$raw_exposure[rows]+
    sinmf_observation_block(sinmf_slice_observation(f$offset,channels),rows)
}

sinmf_stream_operation <- function(task,f,meta,tiny,old_threads,gid) {
  blocks<-split(seq_len(meta$m),ceiling(seq_len(meta$m)/32L));ll<-0
  selected<-task$selected
  if(task$operation=="gradient")gradient<-information<-matrix(0,length(selected),meta$m)
  for(js in blocks) {
    Y<-sinmf_read_matrix_file(file.path(task$directory,"observed.bin"),meta$n,meta$m,channels=js)
    mu<-sinmf_block_mean(f,task$P,task$scale,js)
    W<-sinmf_slice_observation(f$weights,js)
    phi<-sinmf_slice_noise(f$phi,js);theta<-sinmf_slice_noise(f$theta,js)
    if(task$operation=="score")ll<-ll+slcn_loglik(Y,mu,f$kind,phi,theta,as.numeric(W),f$ncores,f$censor)
    if(task$operation=="gradient"&&length(selected)) {
      z<-slcn_pooled_spectra(Y,mu,as(f$H[,selected,drop=FALSE],"dgCMatrix"),
        f$exposure,f$area[selected]*task$scale[selected],f$kind,phi,theta,as.numeric(W),f$ncores,f$censor)
      gradient[,js]<-z$gradient;information[,js]<-z$information
    }
    if(task$operation=="background"&&ncol(f$background)) {
      B<-f$background*f$raw_exposure;A<-f$alpha[,js,drop=FALSE];fixed<-pmax(mu-B%*%A,0)
      precision<-sinmf_weight_matrix(W,nrow(Y))
      if(length(phi)==1L)precision<-precision/phi else precision<-matrix(rep_len(as.numeric(precision),length(Y)),nrow(Y),ncol(Y))/rep(phi,each=nrow(Y))
      z<-if(f$kind==0L)nuisance_gaussian_cpp(B,Y,fixed,precision,A,rep(1L,ncol(B)),f$ncores)else
        if(f$kind%in%1:2)nuisance_poisson_cpp(B,Y,fixed,precision,A,80L,1e-8,f$ncores)else NULL
      if(is.null(z)) {
        z<-slcn_spectral(Y,mu,as(Matrix(B,sparse=TRUE),"dgCMatrix"),t(A),f$kind,phi,theta,
          as.numeric(W),nrow(Y),0L,4L,f$ncores,20L,1e-8,f$censor)
        z<-list(coefficients=t(z$S),background=z$mu-fixed)
      }
      score<-function(mean)slcn_loglik(Y,mean,f$kind,phi,theta,as.numeric(W),f$ncores,f$censor)
      # All channels of a native nuisance update improve independently.
      if(score(fixed+z$background)>=score(mu)-1e-7)f$alpha[,js]<-z$coefficients
    }
  }
  if(task$operation=="score")return(ll)
  if(task$operation=="gradient")return(if(length(selected))list(gradient=gradient,information=information)else NULL)
  f$gid<-gid;f$ncores<-old_threads
  sinmf_checkpoint_sample(f,task$directory,task$operation)
}

# Sum-Poisson proposals are inexpensive; acceptance uses all original ions.
sinmf_aggregate_profile <- function(task,f,meta,tiny,old_threads,gid) {
  Y<-task$aggregate;groups<-task$groups
  P<-aggregateChannels(task$P,groups,meta$mass)
  alpha<-aggregateChannels(f$alpha,groups,meta$mass)
  sizes<-tabulate(sinmf_channel_index(groups,meta$mass),nbins=groups$k)
  O<-if(is.null(dim(f$offset)))outer(rep_len(f$offset,meta$n),sizes)else
    aggregateChannels(f$offset,groups,meta$mass)
  mu<-as.matrix(f$H%*%(P*(f$area*task$scale)))*f$exposure+
    f$background%*%alpha*f$raw_exposure+O
  block<-task$selected
  for(k in block[order(f$area[block]*task$scale[block],decreasing=TRUE)]) {
    at<-which(f$H[,k]>0);if(length(at)<3L)next
    h<-as.numeric(f$H[at,k]);sp<-task$P[k,]*f$area[k]*task$scale[k]
    sp32<-P[k,]*f$area[k]*task$scale[k];W<-sinmf_observation_block(f$weights,at)
    z<-slcn_log_profile(Y[at,,drop=FALSE],mu[at,,drop=FALSE],h,sp32,f$kind,f$phi,f$theta,
      as.numeric(W),length(at),censor=f$censor,maxit=as.integer(task$extra$steps),
      q=as.integer(task$extra$q),nthreads=f$ncores,mode=2L,time=f$aligned_time[at],exposure=f$exposure[at])
    hnew<-z$profile;scale<-sum(hnew*f$quadrature[at]);if(!is.finite(scale)||scale<=0)next
    fullY<-sinmf_read_matrix_file(file.path(task$directory,"observed.bin"),meta$n,meta$m,rows=at)
    oldmu<-sinmf_block_mean(f,task$P,task$scale,seq_len(meta$m),at)
    nextmu<-oldmu+(f$exposure[at]*(hnew-h))%o%sp
    delta<-slcn_loglik(fullY,nextmu,f$kind,f$phi,f$theta,as.numeric(W),f$ncores,f$censor)-
      slcn_loglik(fullY,oldmu,f$kind,f$phi,f$theta,as.numeric(W),f$ncores,f$censor)
    if(is.finite(delta)&&delta>= -1e-7) {
      f$H[at,k]<-hnew/scale;f$area[k]<-f$area[k]*scale
      mu[at,]<-mu[at,,drop=FALSE]+(f$exposure[at]*(hnew-h))%o%sp32
    }
  }
  f$gid<-gid;f$ncores<-old_threads
  sinmf_checkpoint_sample(f,task$directory,"profile")
}
