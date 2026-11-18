test_that("automatic multisample budgets use logical cores and respect worker caps", {
  local_mocked_bindings(detectCores=function(logical=TRUE)if(logical)16L else 8L,
    .package="parallel")
  runtime<-sinmf_runtime(NULL,16L,workers=1L)
  on.exit(sinmf_close_runtime(runtime))
  expect_equal(runtime$budget,list(ncores=16L,workers=1L,threads=16L))
  expect_equal(sinmf_core_count(8L),8L)
  expect_error(sinmf_core_count(0L),"positive integer")
})
