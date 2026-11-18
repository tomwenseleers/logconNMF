test_that("the shared-spectrum sweep uses each sample's full likelihood", {
  n <- 35L
  t <- cumsum(rep(c(.8,1.1,1),length.out=n))
  p <- c(.5,0,.3,.2)
  h <- dnorm(t,mean(t),3)
  h <- h/sum(h)
  w <- outer(seq(.5,1,length.out=n),c(1,2,1,1))
  w[8,1] <- 0
  for(family in c("gaussian","poisson","quasipoisson","nb","nb1","Gamma","Gamma_scale")){
    phi <- switch(family, poisson=1, nb=1, Gamma=.25, 2)
    spec <- slcn_family(family,phi,if(family=="nb")12 else NULL)
    if(is.null(spec$theta))spec$theta <- 1
    desc <- lapply(1:2,function(i){
      area <- 800*i
      off <- outer(seq(.8,1.2,length.out=n),c(2,3,1,2))
      Y <- round(1 + off + (area*h)%o%p)
      if(family %in% c("Gamma","Gamma_scale"))Y[1,1] <- 0
      list(y=Y,time=t,weights=w,offset=off,small_width=5,
           fit=list(H=methods::as(Matrix::Matrix(matrix(h,n,1),sparse=TRUE),"dgCMatrix"),
                    S=matrix(area*p,1),spec=spec,background=matrix(1,n,1),alpha=matrix(1,1,4)))
    })
    fit <- logconNMF(desc,family=family,align=FALSE,precluster=FALSE,ncores=1L,
                     control=list(initial_cycles=0L,cycles=1L,rounds=1L))
    expect_true(fit$batch_accepted, label=family)
    expect_equal(fit$rank,1L)
    expect_equal(fit$L,2L)
    for(path in fit$paths){
      f <- readRDS(path)
      expect_equal(nrow(f$profile_limits),2L)
      z <- as.numeric(f$H[,1]);i <- which(z>0)
      expect_lt(max(diff(diff(log(z[i]))/diff(f$aligned_time[i]))),1e-7)
    }
    expect_true(is.finite(GIC(fit)))
  }
})
