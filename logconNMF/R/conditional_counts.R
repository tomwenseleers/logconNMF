#' @rdname logconNMF
#' @export
conditionalSolverAvailable <- function() is.loaded("_logconNMF_cc_profile_sweep", PACKAGE="logconNMF")

conditional_route <- function(algorithm, family, dispersion, theta, samples) {
  if (algorithm == "classic") return(FALSE)
  unit <- function(x, v) is.null(x) || (!is.character(x) && all(is.finite(x)) && all(x == v))
  supported <- function(s) {
    f <- if(is.character(s$fit)) readRDS(s$fit) else s$fit
    if(!is.null(f$model)) f <- f$model
    w <- if(is.null(s$weights)) f$weights else s$weights
    o <- if(is.null(s$offset)) f$offset else s$offset
    unit(w,1) && (unit(o,0) || (identical(f$algorithm,"conditional") && unit(o,1e-8))) &&
      unit(s$exposure,1) && (is.null(f) || f$spec$kind %in% c(1L,2L))
  }
  valid <- slcn_family(family, dispersion, theta)$kind %in% c(1L, 2L) &&
    all(vapply(samples, supported, logical(1)))
  if (!valid && algorithm == "conditional")
    stop("The conditional solver requires identity-link Poisson or quasipoisson counts, unit weights and exposure, and zero observation offset.")
  valid
}

conditional_windows <- function(h) sinmf_landmarks(h, seq_len(nrow(h)))$windows

conditional_dictionary <- function(X, Y, B, spec, groups, lambda, ncores, control=list()) {
  if (!is.null(groups)) Y <- aggregateChannels(Y, groups)
  X <- as(X, "dgCMatrix")
  B <- as.matrix(B)
  k <- ncol(X)
  alpha_rows <- k + seq_len(ncol(B))
  arguments <- list(X=cbind(X, as(unclass(B), "dgCMatrix")), y=Y,
    family = stats::poisson("identity"), intercept = FALSE, nonnegative = TRUE, grouped = TRUE,
    weights = 1/spec$dispersion, offset = 1e-8, lambda = lambda,
    always.include = alpha_rows, normalize = "none", screen = "none",
    local.search = FALSE, ncores = ncores, group.maxit = 6L)
  settings <- if (is.null(control$initial_solver)) list() else control$initial_solver
  coef <- as.matrix(stats::coef(slcn_dictionary_fit(arguments, settings, !identical(control$initial_fallback,FALSE))))
  s <- coef[seq_len(k), , drop = FALSE]
  live <- which(rowSums(s) > 0 & Matrix::colSums(X) > 0)
  h <- X[, live, drop = FALSE]
  size <- as.numeric(Matrix::colSums(h))
  if (length(size)) h <- h %*% Matrix::Diagonal(x = 1/size)
  list(H = as(h, "dgCMatrix"), S = s[live, , drop = FALSE] * size, alpha = coef[alpha_rows, , drop = FALSE],
       background = B, spec = spec, grouped_dictionary = !is.null(groups))
}

conditional_run <- function(parts, group, price, trace = FALSE, cycles = 12L, exchanges = 2L, ncores = 1L) {
  conditional_native(parts, group, price, trace, cycles, exchanges, sinmf_core_count(ncores))
}

conditional_single <- function(fit, y, time, price, trace) {
  time <- if (is.null(time)) seq_len(nrow(y)) else as.numeric(time)
  if (length(time) != nrow(y) || any(!is.finite(time)) || any(diff(time) <= 0))
    stop("time must contain one increasing finite coordinate per scan.")
  part <- fit[c("H", "S", "alpha", "background", "spec")]
  part$observed <- y
  part$time <- part$aligned <- time
  result <- conditional_run(list(part), seq_len(ncol(y)), price, trace, ncores=fit$ncores)
  f <- result$factors[[1L]]
  counts <- as.numeric(Matrix::colSums(f$H))
  fit$H <- if (length(counts)) as(f$H %*% Matrix::Diagonal(x = 1/counts), "dgCMatrix") else f$H
  fit$S <- result$P * counts
  colnames(fit$S) <- colnames(y)
  fit$alpha <- f$alpha
  fit$offset <- 1e-8
  fit$mu <- as.matrix(fit$H %*% fit$S) + fit$background %*% fit$alpha + fit$offset
  fit$loglik <- slcn_loglik(as.matrix(y), fit$mu, fit$spec$kind, fit$spec$dispersion,
    fit$spec$theta, 1, fit$ncores, fit$control$censor)
  fit$rank <- length(counts)
  fit$conditional_objective <- result$objective
  fit$objective_trace <- result$trace
  fit$ic_df <- result$objective$shape_df + result$objective$spectral_df + fit$nuisance_df
  fit$lambda <- price
  fit$objective <- -2 * fit$loglik + price * fit$ic_df
  fit$algorithm <- "conditional"
  # The input fit's component annotations refer to the start dictionary.
  fit$source_ids <- fit$state_ids <- fit$members <- NULL
  fit <- slcn_temporal_order(fit)
  fit$shape_df <- f$shape_df[order(vapply(seq_along(counts), function(k) which.max(as.numeric(f$H[, k])), 0L))]
  fit
}
