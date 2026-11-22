# Counts are fitted on native scan coordinates. The Python process shares no
# library state with R and the interchange directory belongs to this call.
.conditional_runtime <- new.env(parent = emptyenv())

#' @rdname logconNMF
#' @export
conditionalSolverAvailable <- function() !is.null(conditional_python())

conditional_python <- function() {
  requested <- getOption("logconNMF.python")
  key <- if (is.null(requested)) "auto" else requested
  if (identical(key, .conditional_runtime$key)) return(.conditional_runtime$python)
  conda <- unname(Sys.which("conda"))
  conda_bin <- if (nzchar(conda)) file.path(dirname(dirname(conda)),
    if (.Platform$OS.type == "windows") "python.exe" else "bin/python")
  choices <- if (is.null(requested)) unique(c(conda_bin, unname(Sys.which(c("python3", "python"))))) else requested
  check <- system.file("python", "conditional.py", package = "logconNMF", mustWork = TRUE)
  for (bin in choices) {
    if (!nzchar(bin) || !file.exists(bin) || grepl("Microsoft[/\\\\]WindowsApps", bin, ignore.case = TRUE)) next
    status <- suppressWarnings(tryCatch(system2(bin, c(shQuote(check), "--check"),
      stdout = FALSE, stderr = FALSE), error = function(e) 1L))
    if (identical(status, 0L)) {
      .conditional_runtime$key <- key
      .conditional_runtime$python <- bin
      return(bin)
    }
  }
  NULL
}

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
  if (!valid) return(FALSE)
  if (!is.null(conditional_python())) return(TRUE)
  why <- "The conditional count solver requires Python with numpy, scipy, numba and threadpoolctl. Select it with options(logconNMF.python = '/path/to/python')."
  if (algorithm == "conditional") stop(why, call. = FALSE)
  if (!isTRUE(.conditional_runtime$announced)) {
    message(why, " Using the classic solver.")
    .conditional_runtime$announced <- TRUE
  }
  FALSE
}

conditional_windows <- function(h) {
  result <- matrix(0L, ncol(h), 3L)
  for (j in seq_len(ncol(h))) {
    at <- seq.int(h@p[j] + 1L, h@p[j + 1L])
    v <- h@x[at]
    row <- h@i[at] + 1L
    result[j, ] <- c(range(row[v >= .01 * max(v)]), row[which.max(v)])
  }
  result
}

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

conditional_run <- function(parts, group, price, trace = FALSE, cycles = 12L, exchanges = 2L) {
  tmp <- tempfile("logcon-counts-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  cache <- getOption("logconNMF.cache", tools::R_user_dir("logconNMF", "cache"))
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  cache <- normalizePath(cache, winslash = "/", mustWork = TRUE)
  writeBin(as.integer(group), file.path(tmp, "groups"), endian = "little")
  m <- length(group)
  rows <- lapply(seq_along(parts), function(i) {
    a <- parts[[i]]
    d <- file.path(tmp, i)
    dir.create(d)
    put <- function(name, values) writeBin(values, file.path(d, name), endian = "little")
    y <- a$observed
    if (!is.character(y)) {
      y <- file.path(d, "y")
      sinmf_write_matrix_file(a$observed, y)
    }
    h <- a$H
    put("hp", h@p)
    put("hi", h@i)
    put("hx", h@x)
    put("windows", as.integer(conditional_windows(h)))
    put("B", as.double(a$background))
    put("A", as.double(a$alpha))
    put("S", as.double(a$S))
    put("t", as.double(a$time))
    put("warp", as.double(a$aligned))
    data.frame(n = nrow(h), m = m, b = ncol(a$background), k = ncol(h), q = ncol(a$S),
      phi = a$spec$dispersion, observed = normalizePath(y, winslash = "/", mustWork = TRUE),
      price = price, cycles = cycles, exchanges = exchanges, trace = trace, grouped=isTRUE(a$grouped_dictionary))
  })
  utils::write.csv(do.call(rbind, rows), file.path(tmp, "samples.csv"), row.names = FALSE)
  script <- system.file("python", "conditional.py", package = "logconNMF", mustWork = TRUE)
  log <- file.path(tmp, "solver.log")
  status <- system2(conditional_python(), shQuote(c(script, tmp, cache)), stdout = log, stderr = log)
  if (status != 0L) stop(paste(c("Conditional factorization failed.", readLines(log, warn = FALSE)), collapse = "\n"), call. = FALSE)
  objective <- utils::read.csv(file.path(tmp, "objective.csv"))
  read_values <- function(file, integer = FALSE) {
    size <- if (integer) 4L else 8L
    bytes <- file.size(file)
    if (is.na(bytes) || bytes %% size != 0) stop("Incomplete conditional solver output.")
    readBin(file, if (integer) "integer" else "double", n = bytes/size, size = size, endian = "little")
  }
  K <- objective$rank
  P <- matrix(read_values(file.path(tmp, "spectra")), K, m)
  factors <- lapply(seq_along(parts), function(s) {
    d <- file.path(tmp, s)
    h <- Matrix::sparseMatrix(i = read_values(file.path(d, "fi"), TRUE),
      p = read_values(file.path(d, "fp"), TRUE), x = read_values(file.path(d, "fx")),
      index1 = FALSE, dims = c(nrow(parts[[s]]$H), K))
    list(H = h, alpha = matrix(read_values(file.path(d, "fa")), ncol(parts[[s]]$background), m),
         shape_df = read_values(file.path(d, "df"), TRUE))
  })
  list(P = P, factors = factors, objective = objective,
    trace = if (trace) utils::read.csv(file.path(tmp, "trace.csv")) else NULL)
}

conditional_single <- function(fit, y, time, price, trace) {
  time <- if (is.null(time)) seq_len(nrow(y)) else as.numeric(time)
  if (length(time) != nrow(y) || any(!is.finite(time)) || any(diff(time) <= 0))
    stop("time must contain one increasing finite coordinate per scan.")
  part <- fit[c("H", "S", "alpha", "background", "spec")]
  part$observed <- y
  part$time <- part$aligned <- time
  result <- conditional_run(list(part), seq_len(ncol(y)), price, trace)
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
