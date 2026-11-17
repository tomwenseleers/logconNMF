restart_example <- function() {
  d <- polish_example(2L)
  h <- exp(-(seq_len(nrow(d$f$H))-65)^2/(2*4^2))
  h[c(1:45,86:120)] <- 0
  d$f$H[,2L] <- h/sum(h)
  d$f$S[2L,] <- c(1000,200,10000,5000)
  d$f$offset <- 0;d$f$weights <- 1
  d$f$aggregation <- list(group=c(1L,1L,2L,2L))
  d$Y <- as.matrix(d$f$H%*%d$f$S)+d$f$background%*%d$f$alpha
  d$f <- slcn_sync(d$f,d$Y)
  d
}

test_that("aggregated restart proposals retain the original objective and nuisance terms", {
  d <- restart_example()
  z <- slcn_profile_starts(d$f,d$Y,cycles=2L,knots=32L)
  expect_true(z$profile_starts$aggregated)
  expect_equal(z$profile_starts$screen_channels,2L)
  expect_lte(z$objective,d$f$objective+1e-6)
  expect_equal(dim(z$S),dim(d$f$S))
  expect_identical(z$alpha,d$f$alpha)
  expect_identical(z$offset,d$f$offset)
  expect_identical(z$weights,d$f$weights)
  expect_identical(z$spec,d$f$spec)
  expect_equal(z$mu,as.matrix(z$H%*%z$S)+z$background%*%z$alpha+z$offset,
    tolerance=1e-9)
  expect_equal(as.numeric(Matrix::colSums(z$H)),rep(1,ncol(z$H)),tolerance=1e-9)
  for(k in seq_len(ncol(z$H))) {
    h <- as.numeric(z$H[,k]);at <- which(h>0)
    expect_true(all(diff(at)==1L))
    expect_lte(max(diff(log(h[at]),differences=2L)),1e-6)
  }
  expect_true(all(z$profile_starts$audit$delta[z$profile_starts$audit$accepted]<0))
})

test_that("restart aggregation is restricted to a compatible working model", {
  d <- restart_example();f <- d$f
  expect_equal(slcn_restart_groups(f,d$Y),c(1L,1L,2L,2L))
  for(kind in c(0L,3L,4L,5L,6L)) {
    z <- f;z$spec$kind <- kind
    expect_null(slcn_restart_groups(z,d$Y))
  }
  z <- f;z$offset <- .2
  expect_null(slcn_restart_groups(z,d$Y))
  z <- f;z$weights <- matrix(1,nrow(d$Y),ncol(d$Y))
  expect_null(slcn_restart_groups(z,d$Y))
  z <- f;z$weights <- rep(1,nrow(d$Y))
  expect_null(slcn_restart_groups(z,d$Y))
  z <- f;z$spec$dispersion <- rep(1,ncol(d$Y))
  expect_null(slcn_restart_groups(z,d$Y))
  z <- f;z$aggregation <- NULL
  expect_null(slcn_restart_groups(z,d$Y))
})

# These controlled scores make aggregation look attractive even when the full
# fit is worse. They test the acceptance boundary and cache validity separately
# from the native numerical checks above.
restart_guard <- function(mode=c("reject_entry","reject_control","accept")) {
  mode <- match.arg(mode);d <- restart_example();f <- d$f
  n <- nrow(d$Y);m <- ncol(d$Y);g <- 2L
  fit_score <- function(z)10-sum(as.matrix(z$H)^2)
  f$mu <- matrix(fit_score(f)/m,n,m)
  coarse_calls <- full_calls <- 0L
  local_mocked_bindings(
    slcn_sync=function(f,Y)f,
    slcn_loglik=function(Y,mu,...)mean(mu)*ncol(Y),
    slcn_refine_block=function(z,Y,...) {
      score <- fit_score(z)
      if(ncol(Y)==g)coarse_calls<<-coarse_calls+1L else {
        full_calls<<-full_calls+1L
        if(mode=="reject_entry")score<-score-100
        if(mode=="reject_control"&&full_calls%%2L==0L)score<-score+100
      }
      z$mu <- matrix(score/ncol(Y),nrow(Y),ncol(Y));z
    },.package="logconNMF")
  z <- slcn_profile_starts(f,d$Y,cycles=2L,knots=16L)
  list(fit=z,entry=f,coarse=coarse_calls,full=full_calls)
}

test_that("a favourable aggregated screen cannot bypass either original-data check", {
  for(mode in c("reject_entry","reject_control")) {
    z <- restart_guard(mode)
    expect_true(all(z$fit$profile_starts$audit$start!="none"))
    expect_false(any(z$fit$profile_starts$audit$accepted))
    expect_identical(z$fit$H,z$entry$H)
    expect_gte(z$full,4L)
    expect_equal(z$coarse,7L)
    expect_equal(z$fit$profile_starts$cache_hits,1L)
  }
})

test_that("accepted original-data moves invalidate cached controls", {
  z <- restart_guard("accept")
  expect_true(all(z$fit$profile_starts$audit$accepted))
  expect_equal(z$coarse,8L)
  expect_equal(z$fit$profile_starts$cache_hits,0L)
  expect_equal(z$fit$profile_starts$cache_misses,2L)
})
