test_that("full audits can reactivate a factor and always close the schedule", {
  log <- list()
  step <- function(state, ids, i) { log[[i]] <<- ids; state + 1L }
  schedule <- nmf_sweep_schedule(0L, 1:3, 13L, TRUE, step, identity,
    function(a,b) if (a == 0L) c(1,0,0) else c(0,0,1))
  expect_equal(log[c(1,10,13)], rep(list(1:3),3))
  expect_equal(log[[2]],1L)
  expect_equal(log[[11]],3L)
  expect_equal(schedule$visited,sum(lengths(log)))
  expect_equal(schedule$possible,39L)
  expect_error(nmf_refinement_control(list(refinement_sweeps=1.1)),"integer")
  expect_error(nmf_refinement_control(list(selective=1)),"TRUE or FALSE")
})

test_that("the conditional penalty search returns the nominal objective on both routes", {
  n <- 54L
  time <- seq_len(n)
  X <- convMatrix(dnorm(-5:5,sd=1.5), n)
  H <- as.matrix(X[,c(15,36)])
  Y <- round(2 + H %*% rbind(c(600,300,0,0),c(0,0,500,350)))
  for (family in c("poisson","quasipoisson")) {
    opt <- list(refinement_sweeps=11L, penalty_search=TRUE)
    single <- logconNMF(Y,X,family=family,dispersion=if(family=="poisson") 1 else 2,
      time=time,correlation="none",precluster=FALSE,ncores=1L,control=opt,objective.trace=TRUE)
    batch <- logconNMF(rep(list(list(time=time,y=Y,X=X)),2L),family=family,
      dispersion=if(family=="poisson") 1 else 2,precluster=FALSE,align=FALSE,ncores=1L,
      control=opt,objective.trace=TRUE)
    for (fit in list(single,batch)) {
      expect_equal(fit$refinement$after,fit$conditional_objective$objective,tolerance=1e-8)
      expect_lte(fit$refinement$after,fit$refinement$before + 1e-7)
      expect_lte(fit$refinement$visited,fit$refinement$possible)
      expect_equal(tail(fit$objective_trace$objective,1),fit$refinement$after,tolerance=1e-8)
    }
  }
})

test_that("weighted single and simultaneous fits keep the family likelihood and penalty", {
  n <- 48L
  time <- seq_len(n)
  X <- convMatrix(dnorm(-5:5,sd=1.5),n)
  H <- X[,c(14,31),drop=FALSE]
  H <- H %*% Matrix::Diagonal(x=1/Matrix::colSums(H))
  S <- rbind(c(300,150,0),c(0,100,220))
  offset <- outer(seq(1,2,length.out=n),c(1,2,1))
  weights <- outer(seq(.5,1,length.out=n),c(1,2,1))
  Y <- round(offset + 2 + as.matrix(H%*%S))
  for (family in c("gaussian","poisson","quasipoisson","nb","nb1","Gamma","Gamma_scale")) {
    phi <- switch(family,poisson=1,nb=1,Gamma=.25,2)
    theta <- if (family=="nb") 12
    spec <- slcn_family(family,phi,theta)
    if (is.null(spec$theta)) spec$theta <- 1
    desc <- list(y=Y,time=time,offset=offset,weights=weights,small_width=4,
      fit=list(H=H,S=S,spec=spec,background=matrix(1,n,1),alpha=matrix(2,1,3)))
    opt <- list(refinement_sweeps=2L,penalty_search=TRUE)
    fit <- logconNMF(Y,X,family=family,lambda=3,dispersion=phi,theta=theta,
      weights=weights,offset=offset,background=matrix(1,n,1),ncores=1L,
      precluster=FALSE,control=opt)
    expect_lte(fit$refinement$after,fit$refinement$before+1e-7)
    expect_equal(fit$objective,-2*fit$loglik+fit$lambda*fit$rank)
    joint <- logconNMF(rep(list(desc),2L),family=family,lambda=3,dispersion=phi,theta=theta,
      align=FALSE,ncores=1L,precluster=FALSE,control=c(opt,list(initial_cycles=0L,cycles=1L,rounds=1L)))
    expect_lte(joint$refinement$after,joint$refinement$before+1e-7)
    expect_equal(joint$refinement$after,joint$objective,tolerance=1e-7)
    expect_equal(fitted(joint),fitted(readRDS(file.path(joint$store$path,"model.rds"))))
  }
})

test_that("deleting a weak shared factor retains sample membership and weights", {
  n <- 42L
  X <- convMatrix(dnorm(-4:4), n)
  H <- X[,c(12,30),drop=FALSE]
  S <- rbind(c(200,100),c(0.01,0.02))
  W <- outer(seq(.5,1,length.out=n),c(1,2))
  offset <- matrix(2,n,2)
  spec <- slcn_family("gaussian",2,NULL)
  spec$theta <- 1
  desc <- list(y=as.matrix(H[,1,drop=FALSE]%*%S[1,,drop=FALSE])+3,time=seq_len(n),
    weights=W,offset=offset,small_width=4,
    fit=list(H=H,S=S,spec=spec,background=matrix(1,n,1),alpha=matrix(1,1,2)))
  fit <- logconNMF(rep(list(desc),2),family="gaussian",dispersion=2,lambda=4,
    align=FALSE,ncores=1L,precluster=FALSE,
    control=list(initial_cycles=0L,cycles=0L,rounds=1L,penalty_search=TRUE))
  expect_equal(fit$rank,1L)
  expect_equal(fit$L,2L)
  expect_gt(fit$refinement$before,fit$refinement$after)
  expect_equal(fitted(fit,1L),desc$y,tolerance=1e-4,ignore_attr=TRUE)
})
