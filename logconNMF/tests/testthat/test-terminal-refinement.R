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

test_that("a restricted spectrum sweep holds the background fixed", {
  h <- cbind(c(0,2,3,0,0,0), c(0,0,0,0,0,4))
  p <- rbind(c(.7,.3),c(.1,.9))
  x <- new.env()
  x$h <- h + 0; x$alpha <- matrix(2,1,2); x$z <- matrix(1,6,1)
  x$bg <- x$z %*% x$alpha; x$mu <- x$bg + x$h %*% p
  x$y <- round(2*x$mu); x$t <- x$warp <- as.double(1:6); x$phi <- 1; x$width <- 1
  background <- x$bg + 0; second <- p[2,]
  expect_equal(cc_total_mass(list(x,x), 2:1, 2L), 2*colSums(h)[2:1])
  cc_spectrum_sweep(list(x), p, 1L, 0, FALSE, 2L, background=FALSE)
  expect_equal(x$bg, background); expect_equal(x$alpha, matrix(2,1,2))
  expect_equal(x$h[,2], h[,2]); expect_equal(p[2,], second)
  expect_equal(x$mu, x$bg + x$h %*% p, tolerance = 1e-12)
  cc_spectrum_sweep(list(x), p, 1:2, 0, FALSE, 2L)
  expect_gt(max(abs(x$bg-background)), 0)
  expect_equal(x$mu, x$bg + x$h %*% p, tolerance = 1e-12)
})

test_that("apparent convergence on an early restricted sweep triggers a complete sweep", {
  expect_false(nmf_refinement_control(list())$selective_initial)
  expect_true(nmf_refinement_control(list(penalty_search=TRUE))$selective_initial)
  expect_false(nmf_refinement_control(list(penalty_search=TRUE,refinement_sweeps=100L))$selective_initial)
  expect_true(nmf_refinement_control(list(selective=TRUE))$selective_initial)
  z <- new.env()
  z$data <- list(); z$p <- matrix(1,3,1); z$trace <- FALSE
  z$price <- 1; z$threads <- 1L
  visits <- list()
  local_mocked_bindings(
    cc_D = function(data) 100,
    cc_dimension = function(data,p,D,price) list(objective=D,rank=3L,calls=3L),
    cc_factor_order = function(z,ids) ids,
    cc_refinement_snapshot = function(z) NULL,
    cc_refinement_change = function(before,z) c(1,0,0),
    cc_profile_sweep = function(data,p,order,...) {
      visits[[length(visits)+1L]] <<- order
      if (length(visits) == 1L) 10 else 0
    },
    cc_spectrum_sweep = function(...) 0)
  work <- cc_optimize(z, 8L, "test", selective=TRUE)
  expect_equal(visits,list(1:3,1L,1:3))
  expect_equal(unname(work),c(7,9,3))
})

test_that("screened spectral coordinates agree with independently profiled likelihoods", {
  h <- seq(.1,1,length.out=24)
  offset <- outer(2+h,seq(.5,2,length.out=15))
  y <- round(offset + outer(h,seq(0,1.4,length.out=15)))
  profile <- function(j) {
    score <- function(b) sum(h*(y[,j]/(offset[,j]+b*h)-1))
    if (score(0) <= 0) return(c(coef=0,gain=0))
    b <- uniroot(score,c(0,sum(y[,j])/sum(h)),tol=1e-12)$root
    c(coef=b,gain=2*(sum(y[,j]*log1p(b*h/offset[,j]))-sum(h)*b))
  }
  ref <- vapply(seq_len(ncol(y)),profile,c(coef=0,gain=0))
  for (price in c(0,.12,1,6)) {
    coef <- cc_local_spectrum(y,offset,h,rep(1,ncol(y)),price)
    expect_equal(coef,ifelse(ref["gain",]>price,ref["coef",],0),tolerance=1e-6)
  }
})

test_that("the conditional penalty search returns the nominal objective on both routes", {
  n <- 54L
  time <- seq_len(n)
  X <- convMatrix(dnorm(-5:5,sd=1.5), n)
  H <- as.matrix(X[,c(15,36)])
  Y <- round(2 + H %*% rbind(c(600,300,0,0),c(0,0,500,350)))
  for (family in c("poisson","quasipoisson")) {
    other <- convMatrix(dnorm(-6:6,sd=2),n)
    opt <- list(refinement_sweeps=11L, penalty_search=TRUE, alternative_X=other)
    single <- logconNMF(Y,X,family=family,dispersion=if(family=="poisson") 1 else 2,
      time=time,correlation="none",precluster=FALSE,ncores=1L,control=opt,objective.trace=TRUE)
    batch <- logconNMF(rep(list(list(time=time,y=Y,X=X)),2L),family=family,
      dispersion=if(family=="poisson") 1 else 2,precluster=FALSE,align=FALSE,ncores=1L,
      control=opt,objective.trace=TRUE)
    for (fit in list(single,batch)) {
      expect_equal(fit$refinement$after,fit$restart$before,tolerance=1e-8)
      expect_equal(fit$restart$after,fit$conditional_objective$objective,tolerance=1e-8)
      expect_lte(fit$restart$after,min(fit$restart$before,fit$restart$alternative)+1e-6)
      expect_lte(fit$refinement$after,fit$refinement$before + 1e-7)
      expect_lte(fit$refinement$visited,fit$refinement$possible)
      expect_equal(tail(fit$objective_trace$objective,1),fit$restart$after,tolerance=1e-8)
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
    opt <- list(refinement_sweeps=2L,penalty_search=TRUE,alternative_X=convMatrix(dnorm(-4:4),n))
    fit <- logconNMF(Y,X,family=family,lambda=3,dispersion=phi,theta=theta,
      weights=weights,offset=offset,background=matrix(1,n,1),ncores=1L,
      precluster=FALSE,control=opt)
    expect_lte(fit$refinement$after,fit$refinement$before+1e-7)
    expect_equal(fit$objective,-2*fit$loglik+fit$lambda*fit$rank)
    expect_lte(fit$restart$after,fit$restart$before+1e-7)
    joint <- logconNMF(rep(list(desc),2L),family=family,lambda=3,dispersion=phi,theta=theta,
      align=FALSE,ncores=1L,precluster=FALSE,control=c(opt,list(initial_cycles=0L,cycles=1L,rounds=1L)))
    expect_lte(joint$refinement$after,joint$refinement$before+1e-7)
    expect_equal(joint$restart$after,joint$objective,tolerance=1e-7)
    expect_lte(joint$restart$after,joint$restart$before+1e-7)
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
  desc$y[] <- 3
  desc$fit$S <- desc$fit$S/1e5
  empty <- logconNMF(rep(list(desc),2),family="gaussian",dispersion=2,lambda=4,
    align=FALSE,ncores=1L,precluster=FALSE,
    control=list(initial_cycles=0L,cycles=0L,rounds=0L,penalty_search=TRUE))
  expect_equal(empty$rank,0L)
  expect_equal(dim(empty$areas),c(2L,0L))
  expect_equal(fitted(empty,1L),desc$y,tolerance=1e-5,ignore_attr=TRUE)
})

test_that("an automatic blind count restart can use different sample grids", {
  n <- 90L
  X <- convMatrix(dnorm(-5:5,sd=1.4),n)
  Y <- 2+round(as.matrix(X[,c(14,37,59,75)])%*%rbind(c(800,100),c(0,650),c(60,350),c(220,150)))
  samples <- lapply(c(1,1.01),function(d) list(y=Y,X=X,time=seq_len(n)*d))
  fit <- logconNMF(samples,precluster=FALSE,ncores=1L,dispersion=2,family="quasipoisson",
                   align=FALSE,control=list(alternative_X=TRUE))
  expect_lte(fit$restart$after,min(fit$restart$before,fit$restart$alternative)+1e-6)
  expect_equal(fit$restart$after,fit$conditional_objective$objective,tolerance=1e-8)
  expect_error(logconNMF(samples,control=list(alternative_X=list(X))),"one alternative_X")
})

test_that("native deletion and transfer costs agree with direct original-ion deviances", {
  p <- rbind(c(.1,.7,.2), c(0,.85,.15), c(.7,.1,.2))
  h <- cbind(c(0,1,4,2,0), c(0,3,1,0,0), c(0,0,2,5,1))
  make <- function(h, phi) {
    x <- new.env()
    x$h <- h
    x$phi <- phi
    x$mu <- 2 + h %*% p
    x$y <- round(1 + h %*% p)
    x
  }
  data <- list(make(h, 1), make(cbind(h[,1:2],0), 2.7))
  first <- matrix(c(2L,0L,1L,3L,1L,0L),3L)
  second <- matrix(c(3L,0L,0L,0L,0L,0L),3L)
  fraction <- matrix(0,3L,2L); fraction[1,1] <- .37
  truth <- matrix(NA_real_,3L,2L)
  for (s in 1:2) for (k in 1:3) {
    x <- data[[s]]
    if (!any(x$h[,k] > 0)) next
    spec <- list(rep(0,3))
    j <- first[k,s]; l <- second[k,s]
    if (j) spec[[2]] <- p[j,]
    if (l) spec[[3]] <- (1-fraction[k,s]) * p[j,] + fraction[k,s] * p[l,]
    truth[k,s] <- min(vapply(spec, function(q) {
      m <- x$mu + tcrossprod(x$h[,k], q - p[k,])
      cc_deviance(x$y,m,x$phi) - cc_deviance(x$y,x$mu,x$phi)
    }, 0))
  }
  for (threads in 1:2)
    expect_equal(cc_screen_transfers(data,p,first,second,fraction,threads), truth, tolerance=1e-11)
  expect_error(cc_screen_transfers(data,p,first[-1,,drop=FALSE],second,fraction,1L),"dimensions")
})

test_that("local count search scores a single identity and charges a shared spectrum once", {
  p <- matrix(c(.7,.3),1L)
  h <- matrix(c(0,1,3,1,0),ncol=1)
  make <- function(profile) {
    x <- new.env()
    x$h <- profile; x$t <- x$warp <- as.double(1:5); x$phi <- 1
    x$alpha <- matrix(c(2,2),1,2); x$z <- matrix(1,5,1)
    x$bg <- x$z %*% x$alpha + 1e-8; x$y <- x$bg + 0
    cc_means(x,p); x
  }
  for (N in 1:2) {
    z <- new.env()
    z$p <- p + 0; z$price <- 3; z$trace <- FALSE; z$threads <- 1L
    z$bands <- matrix(c(1,5),1L)
    z$data <- lapply(seq_len(N), function(i) make(h + 0))
    before <- cc_dimension(z$data,z$p,cc_D(z$data),z$price)
    search <- cc_local_search(z,before$objective,TRUE)
    expect_equal(search$removed_by_transfer,N)
    expect_equal(nrow(z$p),0L)
    expect_equal(search$after,0,tolerance=1e-8)
    expect_lt(search$after,search$before)
    for (x in z$data) expect_equal(x$mu,x$bg,tolerance=1e-12)
  }
})

test_that("different sample occurrences of one spectrum can share a single factor", {
  t <- as.double(1:21)
  profile <- round(1e4 * exp(-.3 * (t - 10)^2))
  p <- rbind(c(.6,.4), c(.6,.4))
  z <- new.env()
  z$p <- p + 0; z$price <- 5; z$trace <- FALSE; z$threads <- 1L
  z$bands <- rbind(c(4,16), c(5,17))
  z$data <- lapply(1:2, function(i) {
    x <- new.env()
    x$t <- x$warp <- t; x$phi <- 1
    x$h <- matrix(0,length(t),2); x$h[,i] <- i*profile
    x$z <- matrix(1,length(t),1); x$alpha <- matrix(c(4,3),1,2)
    x$bg <- x$z %*% x$alpha + 1e-8
    cc_means(x,z$p); x$y <- x$mu + 0
    x
  })
  before <- cc_dimension(z$data,z$p,cc_D(z$data),z$price)
  expected_mass <- rowSums(cc_mass(z$data))
  result <- cc_share_nearby_spectra(z)
  cc_trim(z)
  after <- cc_dimension(z$data,z$p,cc_D(z$data),z$price)
  expect_equal(unname(result),c(1,1))
  expect_equal(after$rank,1L); expect_equal(after$calls,2L)
  expect_lt(after$objective,before$objective)
  expect_equal(rowSums(cc_mass(z$data)) / expected_mass,c(1,1),tolerance=.01)
  for (x in z$data) expect_equal(x$mu,x$bg + x$h %*% z$p,tolerance=1e-8)

  # The same RT and a cosine above 0.99 are compatible with two compounds.
  # A sufficiently informative ion must still prevent the proposed merger.
  p <- rbind(c(.7,.3,0), c(.7,.29,.01))
  z$p <- p + 0
  for (i in 1:2) {
    x <- z$data[[i]]
    x$h <- matrix(0,length(t),2); x$h[,i] <- i*profile
    x$alpha <- matrix(2,1,3); x$bg <- x$z %*% x$alpha + 1e-8
    cc_means(x,z$p); x$y <- x$mu + 0
  }
  expect_gt(sum(cc_unit(p)[1,] * cc_unit(p)[2,]), .99)
  expect_equal(unname(cc_share_nearby_spectra(z)),c(1,0))
  expect_equal(z$p,p)
})
