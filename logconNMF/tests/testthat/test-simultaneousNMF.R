test_that("simultaneous fits retain absent occurrences, dimensions and disk-backed reconstruction", {
  set.seed(85);n<-90L;m<-5L;X<-convMatrix(dnorm(-6:6,sd=1.5),n)
  H<-X[,c(25,65)];H<-H%*%Matrix::Diagonal(x=1/colSums(H))
  P<-matrix(runif(2*m),2,m);P<-P/rowSums(P)
  desc<-lapply(1:3,function(s) {
    S<-P*c(16000,12000)*s;if(s==2L)S[2,]<-0
    y<-matrix(rpois(n*m,as.matrix(H%*%S)+2),n,m)
    list(y=y,time=seq_len(n)+s*.2,X=X,background=matrix(1,n,1),small_width=9,
      sample_id=paste0("s",s))
  })
  f<-simultaneousNMF(desc,ncores=1L,control=list(cycles=2L,initial_cycles=1L,rounds=1L,
    initial=list(profile_cycles=0L,original_cleanup=FALSE)))
  on.exit(unlink(f$store$path,recursive=TRUE))
  expect_equal(nrow(f$P),2L);expect_equal(f$L,5L);expect_equal(sum(is.na(f$areas)),1L)
  expect_equal(dim(fitted(f,sample="s2")),c(n,m))
  expect_equal(f$objective,GIC(f))
  expect_equal(fitted(f),fitted(readRDS(file.path(f$store$path,"model.rds"))))
  expect_equal(rowSums(coef(f,normalize="l1")),rep(1,2L),ignore_attr=TRUE)
  expect_true(all(file.exists(f$paths)))
  expect_error(fitted(f,sample="absent"),"Select one sample")
  expect_error(simultaneousNMF(desc,path=f$store$path),"new, empty")
  expect_true(all(f$history$delta[f$history$stage=="simultaneous_cycle"]<=1e-4))
})

test_that("multisample simulations support absence and every signal family", {
  for(family in c("gaussian","poisson","quasipoisson","nb","Gamma","nb1","Gamma_scale")) {
    z<-simulateGCMSSamples(n=c(60L,65L,70L),k=3L,presence=.5,family=family,
      dispersion=if(family=="nb1")2 else 1,seed=8)
    expect_equal(vapply(z$samples,function(x)nrow(x$Y),0L),c(60L,65L,70L))
    expect_true(all(vapply(z$samples,function(x)all(is.finite(x$Y)),FALSE)))
    expect_equal(rowSums(z$areas==0),rowSums(!z$present))
  }
})
