test_that("profile likelihood agrees with direct scalar fits for seven identity families", {
  y <- matrix(c(0, 3, 6, 18, 30), 1)
  p <- c(.05, .1, .15, .3, .4)
  w <- matrix(c(.6, 1, 2, 0, .8), 1)
  off <- matrix(c(.4, .9, 2, 3, 1),1)
  e <- 1.4
  for (family in c("gaussian","poisson","quasipoisson","nb","nb1","Gamma","Gamma_scale")) {
    phi <- switch(family, poisson=1, nb=1, Gamma=.3, 2.4)
    theta <- if(family=="nb") 8 else NULL
    objective <- function(h) {
      mu <- drop(off + e*h*p)
      lp <- switch(family,
        gaussian=dnorm(y, mu, sqrt(phi), log=TRUE),
        poisson=dpois(y, mu, log=TRUE),
        quasipoisson=dpois(y, mu, log=TRUE)/phi,
        nb=dnbinom(y, size=theta, mu=mu, log=TRUE),
        nb1=dnbinom(y, size=mu/(phi-1), mu=mu, log=TRUE),
        {
          shape <- if(family=="Gamma") rep(1/phi,length(p)) else mu/phi
          scale <- mu/shape
          logp <- numeric(length(p))
          zero <- y==0
          logp[zero] <- pgamma(.5,shape[zero],scale=scale[zero],log.p=TRUE)
          logp[!zero] <- dgamma(y[!zero],shape[!zero],scale=scale[!zero],log=TRUE)
          logp
        })
      -sum(w*lp)
    }
    h <- fitLogConcaveProfile(y,p,off,family,exposure=e,weights=w,
                              dispersion=phi,theta=theta,start=2,maxit=40)
    optimum <- optimize(objective,c(.01,120),tol=1e-9)
    expect_equal(as.numeric(h),optimum$minimum,tolerance=3e-4,label=family)
    expect_lte(objective(h)-optimum$objective,1e-7)
  }
  for (family in c("nb1","Gamma_scale")) {
    h <- fitLogConcaveProfile(5,1,offset=5,family=family,dispersion=2,maxit=40)
    best <- optimize(function(mu) if(family=="nb1") -dnbinom(5,size=mu,mu=mu,log=TRUE) else
      -dgamma(5,shape=mu/2,scale=2,log=TRUE),c(5,10))
    expect_equal(h,best$minimum-5,tolerance=.001)
  }
})

test_that("irregular-grid rates obey hard concavity and reconstruct the count model", {
  t <- cumsum(rep(c(.6,1,1.35,.9),length.out=45))
  e <- rep(c(.7,1.2,1.5,.9,1.1),9)
  p <- c(0,.05,.15,.3,.5)
  h <- exp(5-abs(t-23)*ifelse(t<23,.45,.27))
  off <- outer(3+sin(t/5),seq_along(p))
  y <- off+(e*h)%o%p
  for(family in c("gaussian","quasipoisson")) {
    fit <- fitLogConcaveProfile(y,p,off,family,time=t,exposure=e)
    expect_lt(max(abs(fit-h))/max(h),1e-6)
    expect_lt(max(diff(diff(log(fit))/diff(t))),1e-8)
  }
  w <- seq_along(t)/length(t)
  fit <- fitLogConcaveProfile(round(y),p,off,time=t,exposure=e,weights=w)
  scaled <- fitLogConcaveProfile(round(y),3*p,off,time=t,exposure=e,weights=matrix(w,length(t),length(p)))
  expect_equal(3*scaled,fit,tolerance=1e-7)
  expect_equal(fitLogConcaveProfile(c(-2,2,-1),1,family="gaussian"),c(0,2,0),tolerance=1e-6)
  expect_equal(fitLogConcaveProfile(rep(0,5),1,start=rep(1,5)),numeric(5))
  expect_equal(fitLogConcaveProfile(c(3,8),1),c(3,8),tolerance=1e-6)
  expect_error(fitLogConcaveProfile(y,p,family=poisson("log")),"identity")
  expect_error(fitLogConcaveProfile(round(y),p,time=rev(t)),"increasing")
  expect_error(fitLogConcaveProfile(round(y),p,weights=c(1,2)),"dimension|length")
  expect_error(fitLogConcaveProfile(y,p),"counts")
  expect_error(fitLogConcaveProfile(c(1,0.1),1,family="Gamma"),"detection limit")
  expect_error(fitLogConcaveProfile(round(y),p,family="nb"),"theta")
  expect_error(fitLogConcaveProfile(y,p*0),"loading")
})

test_that("all families update noisy irregular profiles on the original scale", {
  set.seed(438)
  n <- 64L; m <- 12L
  p <- seq_len(m)/sum(seq_len(m))
  time <- cumsum(runif(n,.8,1.2))
  e <- runif(n,.6,1.4)
  h <- 1200*exp(-abs(time-32)/ifelse(time<32,4,8))
  old <- 480*exp(-abs(time-26)/5)
  off <- outer(1+sin(time/8)^2,seq_len(m)/4)
  mu <- off+(e*h)%o%p
  for(family in c("gaussian","poisson","quasipoisson","nb","Gamma","nb1","Gamma_scale")){
    phi <- switch(family, gaussian=2,poisson=1,quasipoisson=2,nb=1,Gamma=.07,nb1=2,Gamma_scale=.25)
    y <- switch(family,gaussian=rnorm(n*m,mu,sqrt(phi)),poisson=rpois(n*m,mu),
                quasipoisson=2*rpois(n*m,mu/2),nb=rnbinom(n*m,mu=mu,size=20),
                Gamma=rgamma(n*m,shape=1/phi,scale=mu*phi),nb1=rnbinom(n*m,mu=mu,size=mu/(phi-1)),
                Gamma_scale=rgamma(n*m,shape=mu/phi,scale=phi))
    if(family %in% c("Gamma","Gamma_scale"))y[y<.5]<-0
    fit <- fitLogConcaveProfile(matrix(y,n),p,off,family,time,e,start=old,
             dispersion=phi,theta=if(family=="nb")20 else NULL,maxit=12)
    pos<-which(fit>0)
    expect_true(all(diff(pos)==1));expect_lt(max(diff(diff(log(fit[pos]))/diff(time[pos]))),1e-7)
    expect_gt(sum(h*fit)/sqrt(sum(h*h)*sum(fit*fit)),.96)
  }
  y <- matrix(rpois(128*80,10),128)
  expect_identical(fitLogConcaveProfile(y,rep(1,80),ncores=1),
                   fitLogConcaveProfile(y,rep(1,80),ncores=2))
})
