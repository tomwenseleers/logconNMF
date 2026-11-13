test_that("bundled GC/MS matrices use their nominal mass axes", {
  data("gcmsData1", package = "logconNMF", envir = environment())
  expect_identical(dim(gcmsData1), c(119L, 670L))
  expect_identical(colnames(gcmsData1), as.character(30:699))
  expect_gt(sum(gcmsData1[, 1L]), 0)
  expect_true(is.numeric(attr(gcmsData1, "alkanes_TIC", exact = TRUE)))

  data("gcmsData2", package = "logconNMF", envir = environment())
  expect_identical(dim(gcmsData2), c(1764L, 518L))
  expect_identical(colnames(gcmsData2), as.character(33:550))
})
