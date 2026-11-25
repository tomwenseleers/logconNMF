cc_mass <- function(data) do.call(rbind, lapply(data, function(x) colSums(x$h)))
cc_dimension <- function(data, p, D, price) {
  M <- cc_mass(data)
  live <- colSums(M) > 0
  spectral <- sum(p[live, , drop = FALSE] > 0) - sum(live)
  shape <- sum(vapply(data, function(x) sum(cc_faces(x$t, x$h)), 0))
  data.frame(D = D, objective = D + price * (spectral + shape), rank = sum(live), calls = sum(M > 0),
    spectral_df = spectral, shape_df = shape)
}
cc_D <- function(data) sum(vapply(data, function(x) cc_deviance(x$y, x$mu, x$phi), 0))
cc_history <- function(context, label, cycle, block, D) {
  if (!context$trace)
    return()
  fit <- cc_dimension(context$data, context$p, cc_D(context$data), context$price)
  if (abs(fit$D - D) > 1e-05 + 3e-08 * abs(D))
    stop("Incremental and full original-ion deviances disagree.")
  context$history[[length(context$history) + 1L]] <- cbind(data.frame(stage = label, cycle = cycle,
    block = block, channels = ncol(context$p)), fit)
}
cc_set_h <- function(x, rows, cols, v) cc_write_block(x$h, rows, cols, as.matrix(v))
cc_set_mu <- function(x, rows, cols, m) cc_write_block(x$mu, rows, cols, m)
cc_optimize <- function(z, limit, label, active = FALSE, selective = FALSE) {
  d <- z$data
  p <- z$p
  D <- cc_D(d)
  prior_model <- cc_dimension(d, p, D, z$price)
  cc_history(z, label, 0L, "start", D)
  ids <- use <- seq_len(nrow(p))
  force_full <- TRUE
  visited <- completed <- 0L
  for (it in seq_len(limit)) {
    audit <- !selective || force_full || it %% 10L == 0L || it == limit
    if (audit) use <- ids
    if (audit && selective && it < limit) prior <- cc_refinement_snapshot(z)
    visited <- visited + length(use)
    completed <- it
    full <- length(use) == length(ids)
    D <- D - cc_profile_sweep(d, p, cc_factor_order(z, use), z$bands, z$price, FALSE, z$threads)
    cc_history(z, label, it, "profiles", D)
    D <- D - cc_spectrum_sweep(d, p, cc_factor_order(z, use), z$price, TRUE, z$threads, background = full)
    cc_history(z, label, it, "spectra and background", D)
    support_changed <- FALSE
    if (active) {
      support_gain <- cc_support_sweep(d, p, order(colSums(cc_mass(d))), z$price)
      D <- D - support_gain
      support_changed <- support_gain != 0
      if (it %in% c(2L, 5L)) {
        contraction <- cc_contract(z)
        D <- D - contraction
        support_changed <- support_changed || contraction != 0
      }
      cc_history(z, label, it, "active support", D)
    }
    model <- cc_dimension(d, p, D, z$price)
    support_changed <- support_changed || model$rank != prior_model$rank || model$calls != prior_model$calls
    prev <- prior_model$objective
    now <- model$objective
    if (now > prev + 1e-05 * max(1, prev))
      stop("Count factor update increased the joint objective.")
    small <- prev - now <= 0.001 * prev
    if (small && full) break
    # A support change, even in a full sweep, changes the next conditional
    # problem for neighbouring factors. Recheck them before skipping any.
    force_full <- support_changed || small
    if (audit && selective && it < limit) use <- ids[cc_refinement_change(prior, z) >= 1e-4]
    prior_model <- model
  }
  invisible(c(visited = visited, possible = as.double(completed)*length(ids), sweeps = completed))
}
cc_blocks <- function(z, ids, ions = seq_len(ncol(z$p)), floor_bg = TRUE) {
  blocks <- lapply(seq_along(z$data), function(i)
    cc_local_window(z$data[[i]], z$p, ids, ions, i, floor_bg))
  Filter(Negate(is.null), blocks)
}
cc_local_spectra <- function(bs, vs, ps, ms, price) {
  mass <- Reduce("+", lapply(vs, colSums))
  for (j in order(-mass)) {
    use <- lapply(vs, function(v) which(v[, j] > 0))
    at <- which(lengths(use) > 0)
    if (!length(at)) {
      ps[j, ] <- 0
      next
    }
    os <- ys <- hs <- vector("list", length(at))
    for (l in seq_along(at)) {
      i <- at[l]
      r <- use[[i]]
      b <- bs[[i]]
      v <- vs[[i]][r, j]
      o <- pmax(ms[[i]][r, , drop = FALSE] - tcrossprod(v, ps[j, ]), b$off[r, , drop = FALSE])
      ys[[l]] <- b$y[r, , drop = FALSE]/b$phi
      os[[l]] <- o/b$phi
      hs[[l]] <- v/b$phi
    }
    q <- cc_local_spectrum(do.call(rbind, ys), do.call(rbind, os), unlist(hs, use.names = FALSE),
      ps[j, ], price)
    size <- sum(q)
    for (l in seq_along(at)) {
      i <- at[l]
      r <- use[[i]]
      v <- vs[[i]][r, j]
      ms[[i]][r, ] <- os[[l]] * bs[[i]]$phi + tcrossprod(v, q)
      vs[[i]][r, j] <- v * size
    }
    ps[j, ] <- q/max(size, 1e-100)
  }
  list(v = vs, p = ps, m = ms)
}
cc_coordinate_profiles <- function(b, v, p, m, price, ord = order(-colSums(v)), only_present = FALSE) {
  for (j in ord) {
    if (only_present && !any(v[, j] > 0))
      next
    ion <- which(p[j, ] > 0)
    if (!length(ion))
      next
    q <- p[j, ion]
    o <- pmax(m[, ion, drop = FALSE] - tcrossprod(v[, j], q), b$off[, ion, drop = FALSE])
    v[, j] <- cc_local_profile(b$t, b$y[, ion, drop = FALSE], o, q, v[, j], b$phi, price)
    m[, ion] <- o + tcrossprod(v[, j], q)
  }
  list(v = v, m = m)
}
cc_commit_block <- function(z, bs, ids, fit, ions = seq_len(ncol(z$p))) {
  cc_write_block(z$p, ids, ions, fit$p)
  for (i in seq_along(bs)) {
    b <- bs[[i]]
    x <- z$data[[b$si]]
    cc_set_h(x, b$r, ids, fit$v[[i]])
    cc_set_mu(x, b$r, ions, fit$m[[i]])
  }
}
cc_contract <- function(z) {
  p <- z$p
  band <- z$bands
  mass <- colSums(cc_mass(z$data))
  live <- rowSums(p) > 0
  gain <- 0
  width <- cc_margin(z$data)
  for (k in order(-mass)) {
    if (!live[k])
      next
    c <- which(live & seq_along(live) != k & mass <= mass[k] & band[, 1] <= band[k, 2] + width &
      band[, 2] >= band[k, 1] - width)
    if (!length(c))
      next
    cos <- drop(cc_unit(p[c, , drop = FALSE]) %*% p[k, ])
    for (j in c[utils::head(order(-cos), 2L)]) {
      ids <- c(k, j)
      ion <- which(colSums(p[ids, , drop = FALSE]) > 0)
      bs <- cc_blocks(z, ids, ion, FALSE)
      if (!length(bs))
        next
      weights <- Reduce("+", lapply(bs, function(b) colSums(b$h)/b$phi))
      q <- matrix(drop(weights %*% p[ids, ion, drop = FALSE])/sum(weights), 1L)
      vs <- lapply(bs, function(b) matrix(cc_area_shape(b$t, rowSums(b$h)), ncol = 1L))
      ok <- TRUE
      for (it in 1:2) {
        ms <- vector("list", length(bs))
        for (i in seq_along(bs)) {
          b <- bs[[i]]
          vs[[i]][, 1] <- cc_local_profile(b$t, b$y, b$off, q[1, ], vs[[i]][, 1], b$phi, z$price)
          ms[[i]] <- b$off + vs[[i]] %*% q
        }
        fit <- cc_local_spectra(bs, vs, q, ms, z$price)
        vs <- fit$v
        q <- fit$p
        ms <- fit$m
        if (!any(q > 0)) {
          ok <- FALSE
          break
        }
      }
      if (!ok)
        next
      delta <- cc_trial_cost(bs, vs, q, ms, p[ids, ion, drop = FALSE], z$price)
      if (delta$cost >= -1e-06)
        next
      fit$p <- rbind(q, 0 * q)
      fit$v <- lapply(vs, function(x) cbind(x, 0 * x))
      cc_commit_block(z, bs, ids, fit, ion)
      band[k, ] <- c(min(band[ids, 1]), max(band[ids, 2]))
      live[j] <- FALSE
      mass[k] <- mass[k] + mass[j]
      gain <- gain - delta$d
    }
  }
  z$bands <- band
  gain
}
cc_trial_cost <- function(bs, vs, p, ms, oldp, price) {
  live <- Reduce("+", lapply(vs, colSums)) > 0
  df <- sum(p[live, , drop = FALSE] > 0) - sum(live) - (sum(oldp > 0) - sum(rowSums(oldp) > 0))
  D <- 0
  for (i in seq_along(bs)) {
    b <- bs[[i]]
    D <- D + cc_diff(b$y, ms[[i]], b$old, b$phi)
    df <- df + sum(cc_faces(b$t, vs[[i]])) - sum(cc_faces(b$t, b$h))
  }
  list(cost = D + price * df, d = D)
}
cc_margin <- function(data) {
  widths <- vapply(Filter(function(x) nrow(x$grid) > 0, data), cc_width, 0)
  if (length(widths))
    median(widths) else 0
}
cc_trim <- function(z) {
  keep <- which(colSums(cc_mass(z$data)) > 0 & rowSums(z$p) > 0)
  if (length(keep) != nrow(z$p)) {
    z$p <- z$p[keep, , drop = FALSE]
    z$bands <- z$bands[keep, , drop = FALSE]
    for (x in z$data) x$h <- x$h[, keep, drop = FALSE]
  }
}
