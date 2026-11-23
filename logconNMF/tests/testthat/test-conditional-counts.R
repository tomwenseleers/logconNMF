test_that("conditional factors preserve total ion counts and native-grid shape", {
  expect_true(conditionalSolverAvailable())
  time <- seq(0, 1, length.out=121L)
  X <- convMatrix(convKernel("gaussian", 19L, list(sigma=2)), length(time))
  h <- cbind(exp(-.5*((time-.25)/.035)^2), exp(-.5*((time-.74)/.05)^2))
  p <- rbind(c(.7,.3,0,0), c(0,0,.4,.6))
  y <- round(3 + 2000*h %*% p)
  for (ions in list(1L,1:4)) {
    Y <- y[,ions,drop=FALSE]
    fit <- logconNMF(Y,X,precluster=FALSE,ncores=1L,
      time=time,correlation="none",objective.trace=TRUE)
    expect_identical(fit$algorithm,"conditional")
    expect_gt(fit$rank,0)
    expect_equal(as.numeric(Matrix::colSums(fit$H)),rep(1,fit$rank))
    expect_equal(unname(fitted(fit)),unname(as.matrix(fit$H %*% fit$S)+fit$background%*%fit$alpha+fit$offset))
    expect_lt(sum((fitted(fit)-Y)^2)/sum((Y-mean(Y))^2),.02)
    for(k in seq_len(fit$rank)) {
      v <- as.numeric(fit$H[,k]);i<-which(v>0)
      expect_true(all(diff(i)==1L))
      expect_true(all(diff(diff(log(v[i]))/diff(time[i])) < 1e-5))
    }
    q <- fit$conditional_objective
    expect_equal(q$objective,q$D+fit$lambda*(q$spectral_df+q$shape_df))
    expect_equal(tail(fit$objective_trace$objective,1),q$objective)
  }
})

test_that("shared spectra reconstruct each measured grid and trace the scored model", {
  expect_true(conditionalSolverAvailable())
  set.seed(12)
  n <- 110L
  X <- convMatrix(convKernel("gaussian", 21L, list(sigma=2)), n)
  t <- seq_len(n)
  samples <- lapply(c(0L,2L),function(shift) {
    h<-cbind(2500*exp(-.5*((t-25-shift)/3)^2),1500*exp(-.5*((t-80-shift)/4)^2))
    p<-rbind(c(.5,.35,.15,0),c(0,.1,.15,.75))
    list(y=round(4+h%*%p),time=t*(1+shift/200),X=X)
  })
  for (m in c(1L,4L)) {
    z <- lapply(samples,function(s) {s$y<-s$y[,seq_len(m),drop=FALSE];s})
    ans <- logconNMF(z,precluster=FALSE,ncores=1L,
      family=quasipoisson("identity"),dispersion=2,objective.trace=TRUE)
    expect_identical(ans$algorithm,"conditional")
    expect_equal(rowSums(unname(ans$P)),rep(1,ans$rank))
    expect_equal(ans$L,ans$conditional_objective$calls)
    parallel_fit <- logconNMF(z,precluster=FALSE,ncores=2L,
      control=list(workers=1L),family=quasipoisson("identity"),dispersion=2,
      objective.trace=TRUE)
    expect_equal(parallel_fit$P,ans$P,tolerance=1e-10)
    expect_equal(parallel_fit$areas,ans$areas,tolerance=1e-10)
    expect_equal(parallel_fit$objective_trace,ans$objective_trace,tolerance=1e-8)
    for (i in seq_along(z)) {
      fit <- sampleNMF(ans,i)
      expect_equal(sampleNMF(parallel_fit,i)$H,fit$H,tolerance=1e-10)
      expect_equal(unname(fitted(fit)),unname(fitted(ans,sample=i)),tolerance=1e-8)
      expect_lt(sum((fitted(fit)-z[[i]]$y)^2)/sum(z[[i]]$y^2),.035)
    }
    expect_equal(tail(ans$objective_trace$objective,1),ans$conditional_objective$objective)
  }
  disk <- logconNMF(samples, family=quasipoisson("identity"), dispersion=2,
                    precluster=FALSE, ncores=1L, control=list(cache=FALSE))
  on.exit(unlink(disk$store$path,recursive=TRUE),add=TRUE)
  expect_identical(disk$algorithm,"general")
  expect_null(disk$store$cache)
  expect_true(is.finite(as.numeric(GIC(disk))))
})

test_that("a mixed spectrum returns area to three independent spectra", {
  t <- as.numeric(0:44)
  shape <- dnorm(t,22,4)
  shape <- shape/sum(shape)
  pure <- rbind(c(.7,.3,0,0,0,0),c(0,0,.6,.4,0,0),c(0,0,0,0,.8,.2))
  frac <- c(.25,.35,.4)
  z <- new.env()
  z$p <- rbind(pure,drop(frac%*%pure))
  z$price <- 15
  z$data <- lapply(list(250000*frac,c(80000,0,60000),c(70000,90000,0),c(0,85000,100000)),
    function(area) {
      x<-new.env()
      x$t<-x$warp<-t;x$phi<-2
      h<-tcrossprod(shape,area)
      x$bg<-matrix(7,length(t),ncol(pure));x$y<-x$bg+h%*%pure
      x$h<-cbind(h,0);x
    })
  h<-z$data[[1]]$h
  h[,1:3]<-.12*h[,1:3];h[,4]<-.88*250000*shape;z$data[[1]]$h<-h
  for(x in z$data) cc_means(x,z$p)
  criterion<-function() cc_dimension(z$data,z$p,cc_D(z$data),z$price)$objective
  before<-criterion()
  cc_shared_exchange(z,TRUE);cc_refit_overlaps(z)
  expect_lt(criterion(),before)
  expect_equal(unname(colSums(z$data[[1]]$h)),c(250000*frac,0),tolerance=.005)
  for(x in z$data)expect_true(all(vapply(seq_len(ncol(x$h)),function(j)cc_feasible(x$t,x$h[,j]),TRUE)))
  expect_equal(unname(rowSums(z$p[1:3,,drop=FALSE])),rep(1,3))
  after<-criterion()
  cc_shared_exchange(z,TRUE);cc_refit_overlaps(z)
  expect_equal(criterion(),after,tolerance=1e-8)
})

test_that("native count routing preserves the general identity-link fits", {
  general <- c("gaussian", "gamma", "gamma_scale", "nb1", "nb")
  for (family in general)
    expect_false(conditional_route(family, NULL, NULL, list(list())))
  for (family in c("poisson", "quasipoisson")) {
    expect_true(conditional_route(family, NULL, NULL, list(list())))
    for (s in list(list(weights = c(1, 0)), list(offset = 3), list(exposure = 2)))
      expect_false(conditional_route(family, NULL, NULL, list(s)))
  }
  time <- seq_len(21)
  part <- list(H = as(matrix(0, 21, 0), "dgCMatrix"), S = matrix(0, 0, 2),
    background = matrix(1, 21, 1), alpha = matrix(7, 1, 2), observed = matrix(7, 21, 2),
    time = time, aligned = time, spec = list(dispersion = 2))
  fit <- conditional_run(list(part), 1:2, 15, trace = TRUE, cycles = 2, ncores = 1L)
  expect_equal(fit$objective$rank, 0)
  expect_equal(dim(fit$P), c(0L, 2L))
  expect_equal(fit$factors[[1]]$alpha, part$alpha, tolerance = 1e-7)
})

test_that("local Gram updates satisfy the joint nonnegative least-squares conditions", {
  gram <- Matrix::bandSparse(3L, k=c(-1L, 0L, 1L),
    diagonals=list(c(.4, .2), c(1, 1, 1), c(.4, .2)))
  score <- cbind(c(2, .8, -1), c(-1, 1, 3))
  fit <- cc_gram_nnls(gram, score, matrix(0, 3L, 2L))
  gradient <- as.matrix(gram %*% fit) - score
  expect_true(all(fit >= 0))
  expect_lt(max(abs(gradient[fit > 0])), 1e-7)
  expect_true(all(gradient[fit == 0] >= -1e-7))
})

test_that("local deviance differences match direct count likelihood changes", {
  old <- matrix(c(.5, 1, 2, 4, 8), nrow=1L)
  y <- matrix(c(0, 2, 3, 1, 400), nrow=1L)
  for (ratio in list(c(.1, .5, 1.5, 2, 3), 1 + c(-1, 1, 0, 1, -1) * 1e-12)) {
    mu <- old * rep(ratio, length.out=ncol(old))
    value <- 2 / 1.7 * sum((mu - old) - y * log1p((mu - old) / old))
    expect_equal(cc_diff(y, mu, old, 1.7), value, tolerance=1e-12)
  }
})

test_that("local windows preserve other compounds and background in the offset", {
  x <- new.env()
  x$t <- cumsum(seq(.1,.9,length.out=10))
  x$h <- matrix(0,10,3)
  x$h[3:4,1] <- c(2,4)
  x$h[6:8,2] <- c(5,3,1)
  p <- rbind(c(.7,.3),c(.2,.8),c(0,1))
  x$bg <- matrix(seq(.1,2, length.out=20),10,2)
  x$mu <- x$bg + x$h %*% p
  x$y <- round(x$mu)
  x$phi <- 1.5
  ij <- 3:8
  block <- cc_local_window(x,p,c(2L,1L),2:1,4L,TRUE)
  expect_identical(block$r,ij)
  expect_identical(block$t,x$t[ij])
  expect_equal(block$h,x$h[ij,c(2,1)])
  expect_equal(block$y,x$y[ij,2:1])
  expect_equal(block$off,x$bg[ij,2:1])
  expect_null(cc_local_window(x,p,3L,1L,1L,FALSE))
})

test_that("support deletion checks sample-level shapes and the shared spectrum", {
  make_sample <- function(present) {
    h <- cbind(c(0,3,0,0), c(0,0,3,0))
    z <- new.env()
    z$h <- h
    z$t <- z$warp <- seq_len(4) + 0
    z$phi <- z$width <- 1
    z$bg <- matrix(1,4,2)
    z$mu <- z$bg + h
    z$y <- z$bg
    if (present) z$y[2L,1L] <- 6
    z$z <- matrix(1,4,1)
    z$alpha <- matrix(1,1,2)
    z
  }
  d <- list(make_sample(TRUE),make_sample(FALSE))
  before <- sum(vapply(d,function(x)cc_deviance(x$y,x$mu,x$phi),0))
  p <- diag(2)
  gain <- cc_support_sweep(d,p,1:2,2)
  after <- sum(vapply(d,function(x)cc_deviance(x$y,x$mu,x$phi),0))
  expect_equal(before-after,gain)
  expect_equal(p,diag(c(1,0)))
  expect_equal(unname(colSums(d[[1]]$h)),c(3,0))
  expect_equal(unname(colSums(d[[2]]$h)),c(0,0))
  for (x in d) expect_equal(x$mu,x$bg+x$h%*%p)
})
