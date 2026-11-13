# Internal fitter for sample lists. Public inputs are documented in logconNMF.
simultaneousNMF <- function(samples, family=poisson("identity"), lambda="gic",
                       path=tempfile("simultaneousNMF_"), precluster=32L,
                       dispersion=NULL, theta=NULL, method=c("cosine","hellinger"),
                       align=TRUE, ncores=NULL, control=list(), verbose=FALSE) {
  begun<-proc.time()[["elapsed"]];method<-match.arg(method)
  if(!is.list(samples)||length(samples)<2L)stop("Supply at least two sample descriptions.")
  if(!is.logical(align)||length(align)!=1L||is.na(align))stop("align must be TRUE or FALSE.")
  ncores<-sinmf_core_count(ncores)
  defaults<-list(cycles=10L,initial_cycles=5L,rounds=3L,q=64L,rt_tolerance=NULL,
    alignment_window=NULL,ridge=TRUE,censor=.5,correlation="none",initial=list(),workers=NULL,cache=TRUE)
  if(!is.list(control)||length(control)&&(is.null(names(control))||anyDuplicated(names(control))||
      length(setdiff(names(control),names(defaults)))))stop("Use named simultaneous controls.")
  ctl<-utils::modifyList(defaults,control)
  if(is.list(family)&&!is.null(family$scalar)&&!"censor"%in%names(control))ctl$censor<-family$scalar$censor
  for(nm in c("cycles","initial_cycles","rounds","q")) {
    v<-ctl[[nm]];if(length(v)!=1L||!is.finite(v)||v<(if(nm=="q")3 else 0)||v!=floor(v))stop("Invalid ",nm)
  }
  if(!is.logical(ctl$ridge)||length(ctl$ridge)!=1||is.na(ctl$ridge))stop("Invalid ridge control.")
  if(!is.finite(ctl$censor)||length(ctl$censor)!=1L||ctl$censor<0)stop("Invalid censor.")
  if(dir.exists(path))stop("Use a new, empty store path.")
  dir.create(path,recursive=TRUE);path<-normalizePath(path,winslash="/",mustWork=TRUE)
  if(!is.logical(ctl$cache)||length(ctl$cache)!=1L||is.na(ctl$cache))stop("cache must be TRUE or FALSE.")
  runtime<-sinmf_runtime(ncores,length(samples),ctl$workers)
  on.exit(sinmf_close_runtime(runtime),add=TRUE)
  initial<-system.time(store<-sinmf_prepare(samples,path,family,dispersion,theta,precluster,ncores,ctl,verbose,runtime,lambda))[["elapsed"]]
  state<-sinmf_initialize(store,method,align,ncores,lambda,ctl)
  state$runtime<-runtime;state$parallel<-runtime$budget
  state<-sinmf_optimize(state,cycles=ctl$initial_cycles,q=ctl$q)
  reference<-sinmf_score(state);candidate<-state
  for(round in seq_len(ctl$rounds)) {
    before<-length(unique(candidate$map))
    candidate<-sinmf_batch_spectral_groups(candidate,method,ctl$ridge)
    if(round>1L&&length(unique(candidate$map))==before)break
    candidate<-sinmf_optimize(candidate,cycles=ctl$cycles,q=ctl$q)
    if(verbose)cat("Grouping round",round,":",length(unique(candidate$map)),"compounds; objective",sinmf_score(candidate),"\n")
  }
  candidate$batch_delta<-sinmf_score(candidate)-reference
  accepted<-is.finite(candidate$batch_delta)&&candidate$batch_delta<=0
  if(accepted)state<-candidate
  state$batch_delta<-candidate$batch_delta;state$batch_accepted<-accepted
  state<-sinmf_compact(state)
  state$rank<-nrow(state$P)
  state$loglik<- -sinmf_score(state,penalty=FALSE)/2
  state$ic_df<-(ncol(state$P)-1L)*nrow(state$P)+state$L+state$nuisance_df
  state$objective<- -2*state$loglik+state$lambda*state$ic_df
  state$seconds<-c(initial=initial,simultaneous=proc.time()[["elapsed"]]-begun-initial,total=proc.time()[["elapsed"]]-begun)
  # Only immutable checkpoints of this private store are eligible for cleanup.
  for(d in store$manifest$directory) {
    obsolete<-setdiff(list.files(d,pattern="^fit_.*\\.rds$",full.names=TRUE),state$paths)
    if(length(obsolete))unlink(obsolete)
    meta<-readRDS(file.path(d,"metadata.rds"));meta$initial<-NULL
    saveRDS(meta,file.path(d,"metadata.rds"),compress=FALSE)
  }
  state$runtime<-NULL;state$store$cache<-NULL
  saveRDS(state,file.path(path,"model.rds"))
  state
}

sinmf_edges <- function(time) {
  d<-diff(time);c(time[1L]-d[1L]/2,(utils::head(time,-1L)+utils::tail(time,-1L))/2,utils::tail(time,1L)+utils::tail(d,1L)/2)
}
sinmf_prepare <- function(samples,path,family,dispersion,theta,precluster,ncores,ctl,verbose,runtime=NULL,lambda="gic") {
  groups<-sinmf_common_groups(samples,precluster,ncores)
  first<-samples[[1L]]$y;if(is.character(first))first<-readRDS(first)
  axis<-colnames(first);if(is.null(axis))axis<-as.character(seq_len(ncol(first)))
  tasks<-lapply(seq_along(samples),function(s)list(x=samples[[s]],s=s,path=path,
    family=family,dispersion=dispersion,theta=theta,groups=groups,axis=axis,ctl=ctl,
    threads=if(is.null(runtime))ncores else runtime$budget$threads,verbose=verbose,lambda=lambda))
  results<-sinmf_apply(runtime,tasks,sinmf_prepare_sample)
  manifest<-do.call(rbind,lapply(results,`[[`,"entry"))
  if(anyDuplicated(manifest$sample_id))stop("Use unique sample_id values.")
  saveRDS(manifest,file.path(path,"manifest.rds"),compress=FALSE)
  list(path=path,manifest=manifest,channel_groups=groups,cache=if(ctl$cache)lapply(results,`[[`,"aggregate")else NULL)
}

sinmf_prepare_sample <- function(task) {
    x<-task$x;s<-task$s;path<-task$path;family<-task$family
    dispersion<-task$dispersion;theta<-task$theta;ctl<-task$ctl;precluster<-task$groups
    ncores<-task$threads;verbose<-task$verbose
    if(!is.list(x)||is.null(x$y)||is.null(x$time))stop("Each sample needs y and time.")
    Y<-if(is.character(x$y))readRDS(x$y)else x$y;Y<-as.matrix(Y);storage.mode(Y)<-"double"
    n<-nrow(Y);m<-ncol(Y);time<-as.numeric(x$time)
    if(n<9L||!m||length(time)!=n||any(!is.finite(time))||any(diff(time)<=0)||any(!is.finite(Y)))stop("Invalid sample data/time.")
    mass<-colnames(Y);if(is.null(mass))mass<-as.character(seq_len(m))
    if(!identical(mass,task$axis))stop("Sample channel names and ordering must agree.")
    edges<-if(is.null(x$edges))sinmf_edges(time)else x$edges
    if(length(edges)!=n+1L||any(!is.finite(edges))||any(diff(edges)<=0)||any(time<=utils::head(edges,-1L)|time>=utils::tail(edges,-1L)))stop("Invalid bin edges.")
    E<-if(is.null(x$exposure))rep(1,n)else rep_len(x$exposure,n)
    if(!length(x$exposure)%in%c(0L,1L,n)||any(!is.finite(E)|E<=0))stop("Exposure must be positive and scalar or scan-specific.")
    value<-function(x)if(is.character(x))readRDS(x)else x
    W<-sinmf_observation(value(x$weights),n,m,1);O<-sinmf_observation(value(x$offset),n,m,0)
    if(any(!is.finite(W)|W<0)||!any(W>0)||any(!is.finite(O)|O<0))stop("Invalid weights/offset.")
    B<-if(is.null(x$background))matrix(1,n,1L)else as.matrix(x$background)
    if(nrow(B)!=n||any(!is.finite(B)|B<0))stop("Invalid background basis.")
    f<-if(is.character(x$fit))readRDS(x$fit)else x$fit
    if(!is.null(f$model))f<-f$model
    if(!is.null(f)) {
      if(is.null(x$weights))W<-sinmf_observation(f$weights,n,m,1)
      if(is.null(x$offset))O<-sinmf_observation(f$offset,n,m,0)
      if(any(!is.finite(W)|W<0)||!any(W>0)||any(!is.finite(O)|O<0))stop("Invalid initial-fit weights/offset.")
    }
    if(is.null(f)) {
      X<-x$X
      if(is.null(X))X<-estimatePSF(pmax(Y,0),time=time,family=family,nwindows=1L,precluster=1L,ncores=ncores)$X
      f<-splogconNMF(X,Y,family=family,background=B*E,
        lambda=if(is.numeric(task$lambda))task$lambda*m else task$lambda,precluster=precluster,
        dispersion=if(is.null(x$dispersion))dispersion else x$dispersion,
        theta=if(is.null(x$theta))theta else x$theta,weights=W,offset=O,
        correlation=ctl$correlation,ncores=ncores,control=ctl$initial,verbose=verbose)
    }
    if(nrow(f$H)!=n||ncol(f$S)!=m||ncol(f$H)!=nrow(f$S)||any(!is.finite(f$S)|f$S<0)||
        any(!is.finite(f$H@x)|f$H@x<0))stop("Invalid initial factors.")
    if(is.null(f$spec))stop("Initial fits must contain resolved noise settings in spec.")
    if(f$spec$kind!=0L&&any(Y<0))stop("Nonnegative observations are required by this family.")
    if(f$spec$kind%in%c(1L,3L,5L)&&any(abs(Y-round(Y))>1e-7))stop("Count models require integer observations.")
    if(f$spec$kind%in%c(4L,6L)&&any(Y>0&Y<ctl$censor))stop("Represent values below the Gamma detection limit by zero.")
    id<-if(is.null(x$sample_id))paste0("sample",s)else as.character(x$sample_id)
    if(length(id)!=1||is.na(id)||!nzchar(id))stop("Invalid sample_id.")
    d<-file.path(path,sprintf("sample_%06d",s));dir.create(d)
    sinmf_write_matrix_file(Y,file.path(d,"observed.bin"))
    small<-x$small_width
    if(is.null(small))small<-2*sqrt(2*log(100))*1.6*median(diff(time))
    if(length(small)!=1L||!is.finite(small)||small<=0)stop("small_width must be positive.")
    # Compact initial factors; observed means and reduced-space fits are released.
    f<-f[c("H","S","background","alpha","spec")]
    if(is.null(f$background))f$background<-B*E
    initial<-sinmf_checkpoint_sample(f,d,"initial")
    wfile<-if(length(W)>1L)file.path(d,"weights.rds")else NULL
    ofile<-if(length(O)>1L)file.path(d,"offset.rds")else NULL
    if(!is.null(wfile))saveRDS(W,wfile,compress=FALSE)
    if(!is.null(ofile))saveRDS(O,ofile,compress=FALSE)
    meta<-list(n=n,m=m,time=time,edges=edges,raw_exposure=E,mass=mass,sample_id=id,
      group=if(is.null(x$group))"samples"else as.character(x$group),small_width=small,
      initial=initial,weights=if(is.null(wfile))W else wfile,
      offset=if(is.null(ofile))O else ofile,censor=ctl$censor)
    saveRDS(meta,file.path(d,"metadata.rds"),compress=FALSE)
    aggregate<-if(is.null(precluster))NULL else aggregateChannels(Y,precluster,mass)
    if(!is.null(aggregate)) {
      meta$aggregateChannels<-ncol(aggregate)
      sinmf_write_matrix_file(aggregate,file.path(d,"aggregate.bin"))
      saveRDS(meta,file.path(d,"metadata.rds"),compress=FALSE)
    }
    entry<-data.frame(sample_id=id,directory=d,n=n,m=m,group=meta$group)
    if(verbose)cat("Initialized",id,":",ncol(f$H),"occurrences\n")
    list(entry=entry,aggregate=if(ctl$cache)aggregate else NULL)
}

sinmf_initialize <- function(store,method,align,ncores,lambda,ctl) {
  meta<-lapply(store$manifest$directory,function(d)readRDS(file.path(d,"metadata.rds")))
  # Retain spectra and anchor summaries globally; load profile matrices singly.
  fits<-lapply(meta,function(m) {
    f<-readRDS(m$initial)
    list(S=if(is.null(store$channel_groups))f$S else aggregateChannels(f$S,store$channel_groups,m$mass),spec=f$spec,geometry=sinmf_profile_anchors(f$H,f$S,m$time,m$small_width))
  })
  if(is.null(ctl$rt_tolerance))ctl$rt_tolerance<-max(vapply(meta,`[[`,0,"small_width"))
  if(is.null(ctl$alignment_window))ctl$alignment_window<-5*ctl$rt_tolerance
  for(key in c("rt_tolerance","alignment_window"))if(length(ctl[[key]])!=1L||!is.finite(ctl[[key]])||ctl[[key]]<=0)stop("Invalid ",key)
  alignment<-if(align)sinmf_estimate_alignment(fits,meta,method,lapply(meta,`[[`,"small_width"),ctl$alignment_window,
    geometry=lapply(fits,`[[`,"geometry"))else
    list(maps=lapply(meta,function(x)list(evaluate=identity,method="identity")),reference=1L,method="none")
  N<-sum(store$manifest$n);m<-store$manifest$m[1L]
  factor<-max(vapply(fits,function(f)if(is.null(f$spec$information_factor))1 else f$spec$information_factor,0))
  scalar<-if(is.character(lambda))presetLambda(lambda,p=N,n=N)else lambda
  if(length(scalar)!=1L||!is.finite(scalar)||scalar<0)stop("Invalid lambda.")
  total<-sum(vapply(fits,function(f)nrow(f$S),0L))
  state<-list(store=store,paths=character(length(meta)),P=matrix(0,total,m),map=seq_len(total),
    scale=rep(1,total),serial=0L,lambda=scalar*factor,L=0L,nuisance_df=0L,alignment=alignment,
    history=data.frame(),control=ctl,mass=meta[[1]]$mass,channel_groups=store$channel_groups,n=N,p=N,quasi=FALSE,information_factor=factor)
  for(s in seq_along(meta)) {
    z<-meta[[s]];a<-readRDS(z$initial);fits[s]<-list(NULL);warp<-alignment$maps[[s]]$evaluate
    quadrature<-diff(warp(z$edges));J<-quadrature/diff(z$edges);E<-z$raw_exposure*J
    time<-warp(z$time)
    if(any(!is.finite(E)|E<=0)||any(diff(time)<=0))stop("Alignment must have positive bin widths.")
    H<-a$H;K<-ncol(H);area<-rowSums(a$S);P<-a$S/area
    for(k in seq_len(K)) {
      at<-which(H[,k]>max(H[,k])*1e-10)
      h<-as.numeric(H[at,k])/E[at]
      if(length(at)>2L)h<-sinmf_project_irregular_lc(h,time[at])
      scale<-sum(quadrature[at]*h);H[,k]<-0;H[at,k]<-h/scale;area[k]<-area[k]*scale
    }
    gid<-state$L+seq_len(K);spec<-a$spec
    state$nuisance_df<-state$nuisance_df+as.integer(spec$dispersion_estimated)+as.integer(spec$theta_estimated)
    state$quasi<-state$quasi||spec$kind==2L
    f<-list(H=drop0(as(H,"dgCMatrix")),area=area,gid=gid,aligned_time=time,
      exposure=E,quadrature=quadrature,raw_exposure=z$raw_exposure,
      background=a$background/z$raw_exposure,alpha=a$alpha,offset=z$offset,
      kind=spec$kind,phi=spec$dispersion,theta=spec$theta,spec=spec,censor=z$censor,
      weights=if(is.character(z$weights))readRDS(z$weights)else z$weights,ncores=as.integer(ncores),sample_id=z$sample_id,
      small_width=z$small_width,time=z$time,edges=z$edges)
    if(is.character(f$offset))f$offset<-readRDS(f$offset)
    state$P[gid,]<-P
    state$L<-state$L+K;state<-sinmf_checkpoint(state,s,f,"aligned")
  }
  colnames(state$P)<-state$mass
  state
}

sinmf_sample_index <- function(object,sample) {
  s<-if(is.character(sample))match(sample,object$store$manifest$sample_id)else sample
  if(length(s)!=1L||is.na(s)||s<1||s>length(object$paths)||s!=floor(s))stop("Select one sample by index or name.")
  s
}
#' @noRd
#' @export
fitted.simultaneousNMF <- function(object,sample=1L,...) {
  s<-sinmf_sample_index(object,sample);f<-readRDS(object$paths[s])
  out<-sinmf_mean(f,object);colnames(out)<-object$mass;out
}
#' @noRd
#' @export
coef.simultaneousNMF <- function(object,sample=NULL,normalize=c("none","l1"),...) {
  normalize<-match.arg(normalize)
  if(is.null(sample))return(if(normalize=="l1")object$P else object$consensus_counts)
  f<-readRDS(object$paths[sinmf_sample_index(object,sample)]);ids<-object$map[f$gid]
  out<-object$P[ids,,drop=FALSE]
  if(normalize=="none")out<-out*(f$area*object$scale[ids])
  out
}
#' @noRd
#' @export
logLik.simultaneousNMF <- function(object,...)structure(object$loglik,class="logLik",df=object$ic_df,nobs=object$n,quasi=object$quasi)
#' @noRd
#' @export
nobs.simultaneousNMF <- function(object,...)object$n
#' @noRd
#' @export
AIC.simultaneousNMF <- function(object,...,k=2)-2*object$loglik+k*object$information_factor*object$ic_df
#' @noRd
#' @export
BIC.simultaneousNMF <- function(object,...)AIC(object,k=log(object$n))
#' @noRd
#' @export
GIC.simultaneousNMF <- function(object,...)AIC(object,k=presetLambda("gic",p=object$p,n=object$n))
#' @noRd
#' @export
EBIC.simultaneousNMF <- function(object,...)AIC(object,k=presetLambda("ebic",p=object$p,n=object$n))
#' @noRd
#' @export
MBIC.simultaneousNMF <- function(object,...)AIC(object,k=presetLambda("mbic",p=object$p,n=object$n))
#' @noRd
#' @method print simultaneousNMF
#' @export
print.simultaneousNMF <- function(x,...) {
  cat("Simultaneous log-concave NMF:",length(x$paths),"samples;",nrow(x$P),"shared compounds;",x$L,"occurrences\n")
  cat("Objective:",format(x$objective,digits=9),"; time:",round(x$seconds[["total"]],2),"s\n")
  invisible(x)
}
#' @rdname logconNMF
#' @export
plot.simultaneousNMF <- function(x,type=c("areas","spectra","profiles"),sample=NULL,compounds=seq_len(min(6L,nrow(x$P))),...) {
  type<-match.arg(type)
  if(!length(compounds)||any(!compounds%in%seq_len(nrow(x$P))))stop("Select valid compounds.")
  if(type=="areas")return(invisible(matplot(compounds,t(x$areas[,compounds,drop=FALSE]),type="b",xlab="Compound",ylab="Integrated area",...)))
  axis<-suppressWarnings(as.numeric(x$mass));if(anyNA(axis))axis<-seq_along(x$mass)
  if(type=="spectra")return(invisible(matplot(axis,t(x$P[compounds,,drop=FALSE]),type="h",xlab="Channel",ylab="Spectral fraction",...)))
  samples<-if(is.null(sample))seq_along(x$paths)else sinmf_sample_index(x,sample)
  old<-par(mfrow=c(length(compounds),1L));on.exit(par(old))
  colour<-grDevices::hcl.colors(length(x$paths),"Dark 3")
  for(g in compounds) {
    domain<-range(unlist(lapply(samples,function(s) {
      f<-readRDS(x$paths[s]);k<-which(x$map[f$gid]==g)
      if(!length(k))return(numeric())
      h<-as.numeric(f$H[,k]);f$aligned_time[h>max(h)*1e-4]
    })))
    margin<-max(diff(domain)*.1,.Machine$double.eps)
    plot(domain+c(-margin,margin),c(0,1.05),type="n",xlab="Aligned time",ylab="Relative profile",main=paste("Compound",g),...)
    for(s in samples) {
      f<-readRDS(x$paths[s]);k<-which(x$map[f$gid]==g)
      if(length(k)){h<-as.numeric(f$H[,k]);lines(f$aligned_time,h/max(h),col=colour[s],lwd=1.5)}
    }
    legend("topleft",legend=x$store$manifest$sample_id[samples],col=colour[samples],lty=1,lwd=1.5,bty="n",cex=.8)
  }
  invisible(x)
}

sinmf_noise_pass <- function(state,samples) {
  sinmf_parallel_pass(state,samples,"noise")
}

#' @rdname logconNMF
#' @export
sampleNMF <- function(object,sample=1L) {
  s<-sinmf_sample_index(object,sample);f<-readRDS(object$paths[s]);ids<-object$map[f$gid]
  H<-f$H*f$exposure;sums<-as.numeric(colSums(H))
  S<-object$P[ids,,drop=FALSE]*(f$area*object$scale[ids]*sums)
  if(length(sums))H<-H%*%Diagonal(x=1/sums)
  model<-list(H=as(H,"dgCMatrix"),S=S,background=f$background*f$raw_exposure,alpha=f$alpha,
    offset=f$offset,weights=f$weights,spec=f$spec,n=length(f$time),p=length(f$time),
    ncores=f$ncores,rank=length(ids),compound_id=ids,integrated_area=f$area*object$scale[ids],
    time=f$time,control=list(censor=f$censor),seconds=object$seconds)
  model$loglik<-sinmf_loglik(f,sinmf_sample_data(object$store,s),sinmf_mean(f,object))
  model$nuisance_df<-as.integer(f$spec$dispersion_estimated)+as.integer(f$spec$theta_estimated)
  model$ic_df<-ncol(S)*nrow(S)+model$nuisance_df
  model$lambda<-object$lambda*ncol(S);model$objective<- -2*model$loglik+model$lambda*model$rank
  model<-slcn_temporal_order(model)
  # slcn_temporal_order orders H/S; carry the global IDs through the same ordering.
  apex<-vapply(seq_len(ncol(f$H)),function(k)which.max(f$H[,k]*f$exposure),0L)
  order<-order(apex);model$compound_id<-ids[order];model$integrated_area<-model$integrated_area[order]
  class(model)<-"splogconNMF"
  model
}


# Scalar and row-wise inputs are reused across channels without materialization.
sinmf_observation <- function(x,n,m,default=1) {
  if(is.null(x))return(default)
  if(length(x)==1L)return(as.numeric(x))
  if(length(x)==n&&(is.null(dim(x))||ncol(x)==1L))return(as.numeric(x))
  if(!identical(dim(x),c(n,m)))stop("Observation input needs length n or dimensions n x m.")
  as.matrix(x)
}
sinmf_observation_block <- function(x,rows) {
  if(length(x)==1L)x else if(is.null(dim(x)))x[rows]else x[rows,,drop=FALSE]
}
sinmf_weight_matrix <- function(x,n) {
  if(length(x)==1L)matrix(x,n,1L)else if(is.null(dim(x)))matrix(x,n,1L)else x
}
