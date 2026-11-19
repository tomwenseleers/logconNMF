initial_solver_data <- function() {
  set.seed(472)
  n <- 60L;m <- 3L
  X <- convMatrix(convKernel("gaussian", 9L, list(sigma = 1)), n)
  B <- as.matrix(splines::bs(seq_len(n), df = 4L, intercept = TRUE))
  S <- rbind(c(8000, 4000, 2000), c(1000, 3000, 6000))
  mu <- as.matrix(X[, c(15L, 42L)] %*% S) + 5
  list(X = X, Y = matrix(rpois(n * m, mu), n, m), B = B,
    spec = slcn_family(poisson("identity"), NULL, NULL),
    W = matrix(1, n, m), offset = matrix(.2, n, m), scalar = 3, ncores = 1L)
}

test_that("initial selection can use a matrix-free solver within a zero Gram budget", {
  args <- initial_solver_data()
  constrained <- c(args, list(initial_solver = list(max.memory = 0)))
  expect_error(do.call(slcn_initial, c(constrained, list(initial_fallback = FALSE))),
    "workspace budget")
  expect_message(fallback <- do.call(slcn_initial, constrained),
    "Eigen CG and diagonal preconditioning")
  args$initial_solver <- list(solver = "eigen", method = 3L, precond = 1L, max.memory = 0)
  args$initial_fallback <- FALSE
  explicit <- do.call(slcn_initial, args)
  expect_equal(fallback, explicit, tolerance = 1e-10)
  expect_gt(ncol(explicit$H), 0L)
  expect_true(all(is.finite(explicit$mu)))
  expect_true(all(explicit$S >= 0) && all(explicit$alpha >= 0))
  args$ncores <- 2L
  parallel <- do.call(slcn_initial, args)
  expect_equal(explicit$mu, parallel$mu, tolerance = 1e-8)
})

test_that("public controls reach initialization with and without aggregation", {
  d <- initial_solver_data()
  for (precluster in list(NULL, 2L)) {
    fit <- logconNMF(d$Y, X = d$X, background = d$B, family = quasipoisson("identity"),
      dispersion = 2, offset = d$offset, correlation = "none", precluster = precluster,
      ncores = 1L, control = list(initial_solver = list(method = 3L, precond = 1L,
        max.memory = 0), initial_fallback = FALSE, broad_repeats = 0L,
        full_sweeps = 0L, final_polishes = 0L, profile_cycles = 0L,
        profile_fusion = FALSE, original_cleanup = FALSE))
    expect_s3_class(fit, "logconNMF_fit")
    expect_equal(dim(fitted(fit)), dim(d$Y))
    expect_true(all(is.finite(fitted(fit))))
    expect_equal(fit$control$initial_solver$max.memory, 0)
    if (!is.null(precluster)) expect_equal(fit$aggregation$channels, precluster)
  }
})

test_that("fallback preserves fit settings and propagates unrelated errors", {
  args <- initial_solver_data();calls <- list()
  local_mocked_bindings(l0glmFit = function(...) {
    values <- list(...);calls[[length(calls) + 1L]] <<- values
    if (length(calls) == 1L)
      stop("Poisson Cholesky solve exceeds the Gram/factor workspace budget; use screening.")
    list(coefficients = matrix(1, ncol(values$X), ncol(values$y)))
  })
  args$initial_solver <- list(max.memory = .25, tol = 1e-9)
  expect_message(do.call(slcn_initial, args), "retrying with Eigen CG")
  expect_length(calls, 2L)
  for (name in setdiff(names(calls[[1L]]), c("method", "control.solver")))
    expect_identical(calls[[1L]][[name]], calls[[2L]][[name]])
  expect_identical(calls[[2L]]$method, 3L)
  expect_identical(calls[[2L]]$control.solver$precond, 1L)
  expect_identical(calls[[2L]]$control.solver$max.memory, .25)
  expect_identical(calls[[2L]]$control.solver$tol, 1e-9)
  local_mocked_bindings(l0glmFit = function(...) stop("unrelated numerical error"))
  expect_error(do.call(slcn_initial, args), "unrelated numerical error")
  local_mocked_bindings(l0glmFit = function(...)
    stop("LLT fallback exceeds the Gram/factor workspace budget"))
  args$initial_solver <- list(method = 3L, precond = 1L)
  expect_error(do.call(slcn_initial, args), "LLT fallback exceeds")
})

test_that("simultaneous initialization and refinement retain a zero Gram budget", {
  d <- initial_solver_data()
  samples <- lapply(1:2, function(sample) list(y = d$Y, X = d$X, background = d$B,
    time = seq_len(nrow(d$Y)) + sample / 10, small_width = 9))
  paths <- c(tempfile(), tempfile())
  on.exit(unlink(paths, recursive = TRUE))
  control <- list(workers = 1L, cache = FALSE, initial_cycles = 1L,
    cycles = 1L, rounds = 1L, q = 12L,
    initial = list(initial_solver = list(max.memory = 0)))
  expect_message(fallback <- logconNMF(samples, family = quasipoisson("identity"),
    dispersion = 2, precluster = 2L, correlation = "none", ncores = 1L,
    path = paths[1L], control = control), "Eigen CG and diagonal preconditioning")
  control$initial$initial_solver <- list(method = 3L, precond = 1L, max.memory = 0)
  control$initial$initial_fallback <- FALSE
  control$cache <- TRUE
  explicit <- logconNMF(samples, family = quasipoisson("identity"),
    dispersion = 2, precluster = 2L, correlation = "none", ncores = 1L,
    path = paths[2L], control = control)
  expect_equal(coef(fallback), coef(explicit), tolerance = 1e-8)
  expect_equal(fallback$areas, explicit$areas, tolerance = 1e-8)
  expect_equal(GIC(fallback), GIC(explicit), tolerance = 1e-8)
  expect_null(fallback$store$cache)
  expect_true(any(fallback$history$stage == "simultaneous_cycle"))
  for (sample in 1:2) {
    expect_equal(dim(fitted(fallback, sample = sample)), dim(d$Y))
    expect_true(all(is.finite(fitted(fallback, sample = sample))))
  }
})

test_that("initial solver controls are validated before fitting", {
  expect_error(slcn_initial_solver(list(3L)), "named initial_solver")
  expect_error(slcn_initial_solver(list(unknown = 3L)), "named initial_solver")
  expect_error(slcn_initial_solver(list(solver = "unknown")), "solver")
  expect_error(slcn_initial_solver(list(method = "CG")), "method")
  expect_error(slcn_initial_solver(list(method = 7L)), "method")
  expect_error(slcn_initial_solver(list(method = 3L, precond = 7L)), "precond")
  expect_error(slcn_initial_solver(list(max.memory = -1)), "max.memory")
  d <- initial_solver_data()
  expect_error(splogconNMF(d$X, d$Y, ncores = 1L,
    control = list(initial_fallback = NA)), "initial_fallback")
})
