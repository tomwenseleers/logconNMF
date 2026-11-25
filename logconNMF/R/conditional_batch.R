conditional_workers <- function(samples, ncores) {
  candidates <- Filter(function(s) is.null(s$fit), samples)
  if (!length(candidates)) return(min(ncores,length(samples)))
  # Allow a worker process and workspace for a factorization of the local
  # design. Sparse Gram fill depends on the actual PSF and background.
  columns <- max(vapply(candidates, function(s)
    as.double(if(is.null(s$X)) length(s$time) else ncol(s$X)) + NCOL(s$background), 0))
  per_worker <- 384*1024^2 + 24*columns^2
  available <- cc_available_memory()
  limit <- min(ncores, length(samples))
  if(length(available)==1L && is.finite(available) && available>0)
    limit <- min(limit,max(1,floor(.75*available/per_worker)))
  as.integer(limit)
}

conditional_batch <- function(store, method, align, ncores, lambda, ctl, trace, elapsed, begun, parallel, reference=NULL,
                              refinement=nmf_refinement_control(list())) {
  dirs <- store$manifest$directory
  meta <- lapply(dirs, function(d) readRDS(file.path(d, "metadata.rds")))
  parts <- lapply(seq_along(meta), function(i) {
    x <- meta[[i]]
    a <- readRDS(x$initial)
    a$time <- x$time
    a$observed <- file.path(dirs[i], "observed.bin")
    landmarks <- sinmf_landmarks(a$H, x$time)
    w <- a$windows <- landmarks$windows
    width <- pmax(stats::median(diff(x$time)), x$time[w[, 2L]] - x$time[w[, 1L]])
    a$anchors <- sinmf_profile_anchors(a$H, a$S, x$time, width, landmarks)
    a
  })
  if (is.null(ctl$alignment_window)) {
    width <- if (is.null(ctl$rt_tolerance)) max(vapply(meta, "[[", 0, "small_width")) else ctl$rt_tolerance
    ctl$alignment_window <- 5 * width
  }
  if (length(ctl$alignment_window) != 1L || !is.finite(ctl$alignment_window) || ctl$alignment_window <= 0)
    stop("alignment_window must be positive.")
  # Supplied fits can contain full-channel spectra, whereas PSF regressions can
  # contain grouped spectra. Registration uses the same coordinates for all.
  axis <- meta[[1L]]$mass
  group <- if (is.null(store$channel_groups)) seq_along(axis) else sinmf_channel_index(store$channel_groups, axis)
  geometry <- lapply(parts, "[[", "anchors")
  anchors <- lapply(parts, function(a) {
    p <- a$S
    if (!is.null(store$channel_groups) && ncol(p) == length(axis)) p <- aggregateChannels(p, store$channel_groups, axis)
    list(S = p)
  })
  alignment <- if (align) sinmf_estimate_alignment(anchors, meta, method, lapply(meta, "[[", "small_width"),
    ctl$alignment_window, geometry = geometry) else list(maps = lapply(meta, function(x)
      list(evaluate = identity, method = "identity")), reference = 1L, method = "none")
  for (i in seq_along(parts)) parts[[i]]$aligned <- alignment$maps[[i]]$evaluate(meta[[i]]$time)
  N <- sum(store$manifest$n)
  infl <- max(vapply(parts, function(a) a$spec$information_factor, 0))
  price <- infl * if (is.character(lambda)) L0glm::presetLambda(lambda, p = N, n = N) else lambda
  result <- conditional_run(parts, group, price, trace, ctl$cycles, ncores=ncores, reference=reference,
                            refinement=refinement)
  P <- result$P
  colnames(P) <- axis
  obs <- result$objective
  state <- list(store = store, paths = character(length(parts)), P = P, map = seq_len(nrow(P)),
    scale = rep(1, nrow(P)), serial = 0L, lambda = price, L = obs$calls, rank = nrow(P),
    alignment = alignment, history = if (!is.null(result$trace)) result$trace else data.frame(),
    control = ctl, parallel = parallel, mass = axis, channel_groups = store$channel_groups,
    n = N, p = N, quasi = FALSE, information_factor = infl, algorithm = "conditional",
    conditional_objective = obs, refinement = result$refinement, objective_trace = result$trace, nuisance_df = 0L, loglik = 0)
  for (i in seq_along(parts)) {
    f <- result$factors[[i]]
    meta_i <- meta[[i]]
    spec <- parts[[i]]$spec
    r <- alignment$maps[[i]]$evaluate
    quadrature <- diff(r(meta_i$edges))
    exposure <- quadrature/diff(meta_i$edges)
    live <- which(Matrix::colSums(f$H) > 0)
    h <- f$H[, live, drop = FALSE]
    density <- h/exposure
    area <- as.numeric(Matrix::colSums(density*quadrature))
    if (length(area)) density <- density %*% Matrix::Diagonal(x = 1/area)
    sample <- list(H = as(density, "dgCMatrix"), area = area, gid = live, shape_df = f$shape_df[live],
      aligned_time = parts[[i]]$aligned, exposure = exposure, quadrature = quadrature, raw_exposure = 1,
      background = parts[[i]]$background, alpha = f$alpha, offset = 1e-8, kind = spec$kind,
      phi = spec$dispersion, theta = spec$theta, spec = spec, censor = meta_i$censor, weights = 1,
      ncores = ncores, sample_id = meta_i$sample_id, small_width = meta_i$small_width,
      time = meta_i$time, edges = meta_i$edges)
    state <- sinmf_checkpoint(state, i, sample, "conditional")
    state$loglik <- state$loglik + sinmf_loglik(sample, sinmf_sample_data(store, i), sinmf_mean(sample, state))
    state$nuisance_df <- state$nuisance_df + as.integer(spec$dispersion_estimated) + as.integer(spec$theta_estimated)
    state$quasi <- state$quasi || spec$kind == 2L
  }
  state <- sinmf_compact(state)
  state$ic_df <- obs$shape_df + obs$spectral_df + state$nuisance_df
  state$objective <- -2*state$loglik + price*state$ic_df
  if (!is.null(reference)) {
    state$gic_gap <- cc_gap_traces(result,infl*presetLambda("gic",p=N,n=N))
    state$gic_gap <- gicGap(state)
  }
  total <- proc.time()[["elapsed"]] - begun
  state$seconds <- c(initial = elapsed, simultaneous = total - elapsed, total = total)
  state$store$cache <- NULL
  for (i in seq_along(dirs)) {
    unlink(meta[[i]]$initial)
    meta[[i]]$initial <- NULL
    saveRDS(meta[[i]], file.path(dirs[i], "metadata.rds"), compress = FALSE)
  }
  saveRDS(state, file.path(store$path, "model.rds"))
  state
}
