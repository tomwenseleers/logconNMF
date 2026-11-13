test_that("irregular profile and pooled spectrum updates use the original likelihood", {
  set.seed(88);n<-41L;m<-3L;t<-cumsum(runif(n,.7,1.3));e<-runif(n,.4,1.6)
  h<-dnorm(t,mean(t),4);s<-c(1600,2300,700);mu<-2+e*h%o%s
  Y<-matrix(rpois(length(mu),mu),n,m);W<-matrix(runif(n*m,.5,1.5),n,m)
  for(kind in 0:6) {
    phi<-if(kind==4L).3 else 3
    run<-function(nt)slcn_log_profile(Y,mu,h,s,kind,phi,10,as.numeric(W),n,
      q=24L,nthreads=nt,time=t,exposure=e)
    a<-run(1L)
    expect_identical(a,run(2L))
    expect_lt(max(diff(diff(log(a$profile))/diff(t))),1e-6)
    expect_true(all(a$history<=1e-7))
    z<-slcn_pooled_terms(Y,mu,e*h,kind,phi,10,W)
    expect_identical(z,slcn_pooled_terms(Y,mu,e*h,kind,phi,10,W,2L))
    eps<-1e-3;d<-e*h%o%c(0,eps,0)
    ll<-function(v)slcn_loglik(Y,v,kind,phi,10,as.numeric(W),1L,.5)
    expect_lt(abs(z$gradient[2]+(ll(mu+d)-ll(mu-d))/(2*eps)),1e-5)
  }
  a<-slcn_log_profile(Y,mu,h,s,1L,1,1,1,n,q=24L)
  b<-slcn_log_profile(Y,mu,h,s,1L,1,1,1,n,q=24L,time=seq_len(n),exposure=rep(1,n))
  expect_identical(a,b)
})

test_that("batched spectral scores equal profile-wise scores with shared rows", {
  set.seed(7);n<-55L;m<-4L
  H<-Matrix::Matrix(outer(seq_len(n),c(15,40),function(i,j)exp(-((i-j)/4)^2)),sparse=TRUE)
  H@x[H@x<1e-12]<-0;H<-as(Matrix::drop0(H),"dgCMatrix")
  e<-runif(n,.8,1.2);a<-c(200,300);S<-matrix(runif(2*m),2,m)
  mu<-as.matrix((H*e)%*%(S*a))+2;Y<-matrix(rpois(n*m,mu),n,m);w<-runif(n,.7,1.3)
  for(kind in 0:6) {
    z<-slcn_pooled_spectra(Y,mu,H,e,a,kind,.5,5,w,1L)
    expect_identical(z,slcn_pooled_spectra(Y,mu,H,e,a,kind,.5,5,w,2L))
    for(k in 1:2) {
      old<-slcn_pooled_terms(Y,mu,as.numeric(H[,k])*e*a[k],kind,.5,5,matrix(w,n,1L),1L)
      expect_equal(z$gradient[k,],old$gradient,tolerance=1e-13)
      expect_equal(z$information[k,],old$information,tolerance=1e-13)
    }
    expect_identical(slcn_loglik(Y,mu,kind,.5,5,w,1L,.5),
      slcn_loglik(Y,mu,kind,.5,5,as.numeric(matrix(w,n,m)),1L,.5))
  }
})
