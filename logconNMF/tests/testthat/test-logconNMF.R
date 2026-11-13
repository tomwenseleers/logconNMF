test_that("one matrix uses the single-sample penalty and extractors", {
  set.seed(86)
  X <- convMatrix(dnorm(-5:5, sd = 1.3), 60L)
  S <- matrix(c(2600, 900, 1100, 2200, 480, 1700), 2L, 3L)
  mu <- as.matrix(X[, c(18L, 43L)] %*% S) + 1
  Y <- matrix(rpois(length(mu), mu), nrow(mu), dimnames = list(NULL, c("43", "57", "71")))
  B <- matrix(1, nrow(Y), 1L)
  direct <- splogconNMF(X, Y, lambda = 2 * ncol(Y), background = B, ncores = 1L)
  fit <- logconNMF(Y, X = X, lambda = 2, background = B, ncores = 1L)
  preset <- logconNMF(Y, X = X, lambda = "aic", background = B, ncores = 1L)
  expect_s3_class(fit, "logconNMF_fit")
  expect_identical(fit$H, direct$H)
  expect_identical(coef(fit), coef(direct))
  expect_identical(coef(fit), coef(preset))
  expect_identical(fitted(fit, sample = "sample1"), fitted(direct))
  expect_identical(GIC(fit), GIC(direct))
  expect_equal(unname(rowSums(coef(fit, normalize = "l1"))), rep(1, fit$rank))
  expect_error(coef(fit, sample = 2L), "Select sample 1")
  empty <- logconNMF(Y * 0, X = X, ncores = 1L, correlation = "none")
  expect_equal(dim(coef(empty, normalize = "l1")), c(0L, ncol(Y)))
})

test_that("simultaneous calls retain sample names, file-backed inputs and joint scores", {
  set.seed(87)
  X <- convMatrix(dnorm(-5:5, sd = 1.3), 60L)
  S <- matrix(c(2300, 800, 1100, 2100, 480, 1600), 2L, 3L)
  mu <- as.matrix(X[, c(18L, 43L)] %*% S) + 1
  y <- matrix(rpois(length(mu), mu), nrow(mu))
  ypath <- tempfile(fileext = ".rds")
  saveRDS(y, ypath)
  samples <- list(first = list(y = y, X = X, time = seq_len(nrow(y))),
    second = list(y = ypath, X = X, time = seq_len(nrow(y)) + .1))
  for (name in names(samples)) samples[[name]]$small_width <- 8
  paths <- c(tempfile(), tempfile())
  on.exit(unlink(c(paths, ypath), recursive = TRUE))
  ctl <- list(cycles = 1L, initial_cycles = 1L, rounds = 1L,
    initial = list(profile_cycles = 0L, original_cleanup = FALSE))
  legacy <- samples
  for (name in names(legacy)) legacy[[name]]$sample_id <- name
  reference <- simultaneousNMF(legacy, path = paths[1], ncores = 1L, lambda = 2, control = ctl)
  fit <- logconNMF(samples, path = paths[2], ncores = 1L, lambda = 2, control = ctl)
  expect_s3_class(fit, "logconNMF_fit")
  expect_equal(fit$rank, nrow(reference$P))
  expect_identical(fit$areas, reference$areas)
  expect_identical(coef(fit, normalize = "l1"), coef(reference, normalize = "l1"))
  expect_identical(fitted(fit, sample = "second"), fitted(reference, sample = 2L))
  expect_identical(GIC(fit), GIC(reference))
  restored <- readRDS(file.path(fit$store$path, "model.rds"))
  expect_equal(restored$rank, fit$rank)
  expect_identical(fitted(restored, sample = "second"), fitted(fit, sample = "second"))
})

test_that("arguments for the wrong fitting path give explicit errors", {
  y <- matrix(0, 20L, 2L)
  sample <- list(y = y, time = seq_len(nrow(y)))
  expect_error(logconNMF(y), "Supply X")
  expect_error(logconNMF(y, X = y, path = "unused"), "apply to simultaneous fits")
  expect_error(logconNMF(list(sample)), "at least two")
  expect_error(logconNMF(list(sample, sample), offset = 1), "in the sample descriptions")
  expect_error(logconNMF(list(sample, list(y = y))), "needs y and time")
  expect_error(logconNMF(y, lambda = c(1, 2)), "one finite")
})
