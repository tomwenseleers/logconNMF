test_that("cached count scores retain absolute constants, weights and boundary values", {
  set.seed(181119);n<-17L;m<-4L
  mu<-matrix(runif(n*m,.1,1000),n,m);Y<-matrix(as.double(rpois(n*m,mu)),n,m)
  mu[1,1]<-Y[1,1]<-0;Y[2,1]<-2
  for(kind in 1:2)for(W in list(1,runif(n),matrix(runif(n*m),n,m)))for(nt in c(1L,2L)) {
    phi<-seq(1,2,length.out=m)
    plain<-slcn_loglik(Y,mu,kind,phi,1,as.numeric(W),nt,.5)
    cached<-Y;attr(cached,"sinmf_log_factorial")<-slcn_log_factorial(Y,nt)
    expect_identical(slcn_loglik(cached,mu,kind,phi,1,as.numeric(W),nt,.5),plain)
    zero<-mu;zero[2,1]<-0
    expect_identical(slcn_loglik(cached,zero,kind,phi,1,as.numeric(W),nt,.5),
      slcn_loglik(Y,zero,kind,phi,1,as.numeric(W),nt,.5))
  }
  bad<-Y;attr(bad,"sinmf_log_factorial")<-matrix(0,1,1)
  expect_error(slcn_loglik(bad,mu,1L,1,1,1),"normalizer dimensions")
  expect_error(slcn_log_factorial(matrix(-1,2,2)),"nonnegative")
  expect_error(slcn_loglik(Y,mu[-1,,drop=FALSE],1L,1,1,1),"dimensions")
})

test_that("one bounded count cache retains slices and is released on exit", {
  cache<-get(".sinmf_cache",asNamespace("logconNMF"));previous<-as.list(cache)
  runtime<-sinmf_runtime(1L,1L,cache=TRUE)
  on.exit(if(!is.null(runtime))sinmf_close_runtime(runtime),add=TRUE)
  directory<-tempfile("coupled-cache-");dir.create(directory)
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  A<-matrix(as.double(seq_len(60)),20,3);B<-A+7
  paths<-file.path(directory,c("a","b"),"observed.bin")
  for(d in dirname(paths))dir.create(d)
  sinmf_write_matrix_file(A,paths[1]);sinmf_write_matrix_file(B,paths[2])
  sinmf_cache_sample(2L)
  for(rows in list(1:20,c(19L,3L,5L)))for(channels in list(1:3,c(3L,1L))) {
    z<-sinmf_read_matrix_file(paths[1],20L,3L,channels,rows)
    expect_identical(as.numeric(z),as.numeric(A[rows,channels,drop=FALSE]))
    expect_identical(attr(z,"sinmf_log_factorial"),
      slcn_log_factorial(A)[rows,channels,drop=FALSE])
  }
  z<-sinmf_read_matrix_file(paths[2],20L,3L)
  expect_identical(as.numeric(z),as.numeric(B))
  expect_true(grepl(paths[2],cache$key,fixed=TRUE))
  expect_equal(length(cache$Y),length(B))
  expect_null(sinmf_cached_read(paths[1],10000000L,2L,1L,1L))
  expect_error(sinmf_read_matrix_file(paths[1],20L,3L,4L),"Invalid matrix slice")
  sinmf_cache_sample(0L)
  expect_null(cache$Y)
  sinmf_write_matrix_file(-A,paths[1])
  z<-sinmf_read_matrix_file(paths[1],20L,3L)
  expect_identical(z,-A);expect_null(attr(z,"sinmf_log_factorial"))
  sinmf_close_runtime(runtime);runtime<-NULL
  expect_identical(as.list(cache),previous)
  expect_error((function() {
    run<-sinmf_runtime(1L,1L,cache=TRUE);on.exit(sinmf_close_runtime(run))
    sinmf_cache_sample(1L);sinmf_read_matrix_file(paths[2],20L,3L)
    stop("forced runtime exit")
  })(),"forced runtime exit")
  expect_identical(as.list(cache),previous)
})
