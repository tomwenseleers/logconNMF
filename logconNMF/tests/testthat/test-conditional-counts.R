test_that("conditional factors preserve total ion counts and native-grid shape", {
  skip_if_not(conditionalSolverAvailable())
  time <- seq(0, 1, length.out=121L)
  X <- convMatrix(convKernel("gaussian", 19L, list(sigma=2)), length(time))
  h <- cbind(exp(-.5*((time-.25)/.035)^2), exp(-.5*((time-.74)/.05)^2))
  p <- rbind(c(.7,.3,0,0), c(0,0,.4,.6))
  y <- round(3 + 2000*h %*% p)
  for (ions in list(1L,1:4)) {
    Y <- y[,ions,drop=FALSE]
    fit <- logconNMF(Y,X,precluster=FALSE,ncores=1L,algorithm="conditional",
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
  skip_if_not(conditionalSolverAvailable())
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
    ans <- logconNMF(z,algorithm="conditional",precluster=FALSE,ncores=1L,
      family=quasipoisson("identity"),dispersion=2,objective.trace=TRUE)
    expect_identical(ans$algorithm,"conditional")
    expect_equal(rowSums(unname(ans$P)),rep(1,ans$rank))
    expect_equal(ans$L,ans$conditional_objective$calls)
    parallel_fit <- logconNMF(z,algorithm="conditional",precluster=FALSE,ncores=2L,
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
  expect_error(logconNMF(samples,algorithm="conditional",family=gaussian()),"requires")
})

test_that("a mixed spectrum can return area to three independent spectra", {
  skip_if_not(conditionalSolverAvailable())
  script <- test_path("fixtures", "mixed-spectrum.py")
  modules <- system.file("python", package="logconNMF", mustWork=TRUE)
  log <- system2(conditional_python(), shQuote(c(script, modules, tempdir())),
                 stdout=TRUE, stderr=TRUE)
  expect_true(is.null(attr(log, "status")), info=paste(log, collapse="\n"))
})
