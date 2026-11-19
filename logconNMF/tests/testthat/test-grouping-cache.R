ridge_full_grid <- function(profiles,times,areas,tolerance=.7,max_shift=.12) {
  o<-order(areas);coverage<-raw_coverage<-shift<-numeric(length(o)-1L)
  for(j in seq_along(coverage)) {
    a<-o[j];b<-o[j+1L];grid<-sort(unique(c(times[[a]],times[[b]])))
    low<-approx(times[[a]],profiles[[a]]*areas[a],grid,yleft=0,yright=0)$y
    high<-approx(times[[b]],profiles[[b]]*areas[b],grid,yleft=0,yright=0)$y
    integral<-function(v)sum(diff(grid)*(head(v,-1L)+tail(v,-1L))/2)
    raw_coverage[j]<-integral(pmin(low,high))/max(1e-20,integral(low))
    start_a<-times[[a]][which(profiles[[a]]>max(profiles[[a]])*.01)[1L]]
    start_b<-times[[b]][which(profiles[[b]]>max(profiles[[b]])*.01)[1L]]
    shift[j]<-max(-max_shift,min(max_shift,start_a-start_b))
    moved<-approx(times[[b]]+shift[j],profiles[[b]]*areas[b],grid,yleft=0,yright=0)$y
    coverage[j]<-max(raw_coverage[j],integral(pmin(low,moved))/max(1e-20,integral(low)))
  }
  list(pass=all(coverage>=tolerance),coverage=coverage,raw_coverage=raw_coverage,shift=shift)
}

test_that("local ridge integrals retain full-grid boundary and irregular-grid scores", {
  set.seed(182611)
  for(i in 1:80) {
    n<-sample(50:200,1);m<-sample(50:200,1)
    ta<-cumsum(runif(n,.008,.012));tb<-cumsum(runif(m,.008,.012))
    if(i%%5==0)tb<-tb+max(ta)
    ca<-sample(seq_along(ta),1);cb<-sample(seq_along(tb),1)
    ha<-dnorm((ta-ta[ca])/.06);hb<-dnorm((tb-tb[cb])/.09)
    ha[abs(ta-ta[ca])>.3]<-0;hb[abs(tb-tb[cb])>.4]<-0
    profiles<-list(ha/sum(ha),hb/sum(hb));times<-list(ta,tb);areas<-exp(rnorm(2))
    expected<-ridge_full_grid(profiles,times,areas)
    actual<-sinmf_ridge_compatibility(profiles,times,areas)
    expect_identical(actual$pass,expected$pass)
    expect_equal(actual,expected,tolerance=1e-12)
  }
})

test_that("cached grouping preserves rejected-edge ordering and sparse factors", {
  store<-tempfile();on.exit(unlink(store,recursive=TRUE))
set.seed(12618);n <- 120L;t <- seq_len(n)/60
P <- rbind(c(.8,.18,.02),c(.79,.19,.02),c(.8,.18,.02))
area <- c(10000,30000,30000);sigma <- c(3,5,.5);center <- c(50,68,50)
descriptions <- lapply(1:3,function(s) {
  h <- dnorm((seq_len(n)-center[s])/sigma[s]);h[abs(seq_len(n)-center[s])>10*sigma[s]] <- 0;h <- h/sum(h)
  f <- list(H=Matrix(matrix(h,n,1L),sparse=TRUE),S=matrix(P[s,]*area[s],1L),
    background=matrix(1,n,1L),alpha=matrix(1,1L,3L),spec=slcn_family(poisson('identity')))
  Y <- matrix(rpois(n*3,h%o%(P[s,]*area[s])+1),n,3,dimnames=list(NULL,as.character(40:42)))
  f$spec <- slcn_noise(Y,slcn_family(poisson('identity')),1,1L,'none',.5)$spec
  list(y=Y,fit=f,time=t,sample_id=paste0('s',s),small_width=.1)
})
state <- logconNMF(descriptions,align=FALSE,path=store,
  precluster=FALSE,ncores=1L,control=list(cycles=0L,initial_cycles=0L,rounds=0L,rt_tolerance=.5))

  expected_P<-state$P;expected_map<-state$map
  z<-sinmf_batch_spectral_groups(state)
  expect_identical(state$P,expected_P);expect_identical(state$map,expected_map)
  expect_equal(z$batch_grouping[c("a","b","accepted")],data.frame(
    a=c(1L,2L,1L,2L),b=c(2L,3L,3L,1L),accepted=c(FALSE,FALSE,TRUE,FALSE)))
  expect_equal(length(unique(z$map)),2L)
  expect_gt(z$grouping_cache_bytes,0)
  expect_lt(z$grouping_cache_bytes,2*1024^2)
  expect_null(z$geometry)
})


test_that("original profile buffers are capped per worker without large allocation", {
  calls<-0L;expected<-matrix(1:4,2L)
  previous<-as.list(.sinmf_cache)
  on.exit(sinmf_cache_restore(previous))
  local_mocked_bindings(sinmf_read_matrix_file=function(path,n,m) {
    calls<<-calls+1L
    expect_equal(basename(path),"observed.bin")
    expected
  })
  sinmf_cache_reset(FALSE)
  expect_null(sinmf_original_cache(tempdir(),10000L,550L))
  expect_equal(calls,0L)
  sinmf_cache_reset(TRUE)
  expect_identical(sinmf_original_cache(tempdir(),10000L,550L),expected)
  expect_identical(sinmf_original_cache(tempdir(),2^20,8L),expected)
  expect_null(sinmf_original_cache(tempdir(),2^20+1,8L))
  expect_null(sinmf_original_cache(tempdir(),1e6,1e6))
  expect_equal(calls,2L)
})
