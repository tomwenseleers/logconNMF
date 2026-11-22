cc_choices <- function(k, c, unit) {
  if (!length(c))
    return(list())
  cos <- drop(unit[c, , drop = FALSE] %*% unit[k, ])
  top <- utils::head(order(-cos), 8L)
  out <- list()
  if (cos[top[1]] > 0.9)
    out <- list(c[top[1]])
  best <- 0.9^2
  pair <- NULL
  if (length(top) > 1L)
    for (a in seq_len(length(top) - 1L)) for (b in seq.int(a + 1L, length(top))) {
      i <- top[a]
      j <- top[b]
      cross <- sum(unit[c[i], ] * unit[c[j], ])
      det <- 1 - cross^2
      if (det <= 1e-10)
        next
      w <- c(cos[i] - cross * cos[j], cos[j] - cross * cos[i])/det
      score <- sum(w * cos[c(i, j)])
      if (min(w) > 0 && score > best) {
        best <- score
        pair <- c[c(i, j)]
      }
    }
  if (length(pair))
    out[[length(out) + 1L]] <- pair
  out
}
cc_exchange_trial <- function(z, si, k, rec) {
  x <- z$data[[si]]
  ids <- c(k, rec)
  ions <- which(colSums(z$p[ids, , drop = FALSE]) > 0)
  pos <- which(rowSums(x$h[, ids, drop = FALSE]) > 0)
  r <- seq.int(pos[1], utils::tail(pos, 1))
  old <- x$mu[r, ions, drop = FALSE]
  prior <- x$h[r, ids, drop = FALSE]
  off <- pmax(old - prior %*% z$p[ids, ions, drop = FALSE], 1e-08)
  b <- list(t = x$t[r], y = x$y[r, ions, drop = FALSE], off = off, phi = x$phi)
  p <- z$p[rec, ions, drop = FALSE]
  v <- prior[, -1, drop = FALSE] + prior[, 1]/length(rec)
  v <- cc_scan_unmix(b$y, off, p, v, matrix(TRUE, nrow(v), ncol(v)), x$phi)
  for (j in seq_along(rec)) v[, j] <- cc_area_shape(b$t, v[, j])
  upd <- cc_coordinate_profiles(b, v, p, off + v %*% p, z$price)
  dd <- cc_diff(b$y, upd$m, old, x$phi)
  list(cost = dd + z$price * (sum(cc_faces(b$t, upd$v)) - sum(cc_faces(b$t, prior))), d = dd, v = upd$v,
    m = upd$m, r = r, ion = ions, rec = rec)
}
cc_exchange <- function(z) {
  p <- z$p
  m <- cc_mass(z$data)
  unit <- cc_unit(p)
  margin <- cc_margin(z$data)
  n <- 0L
  for (k in order(colSums(m > 0), colSums(m))) {
    members <- which(m[, k] > 0)
    if (!length(members))
      next
    band <- z$bands
    c <- which(colSums(m) > 0 & band[, 1] <= band[k, 2] + margin & band[, 2] >= band[k, 1] - margin)
    choices <- cc_choices(k, setdiff(c, k), unit)
    if (!length(choices))
      next
    offers <- lapply(members, function(si) {
      fits <- lapply(choices, function(u) cc_exchange_trial(z, si, k, u))
      fits[[which.min(vapply(fits, "[[", 0, "cost"))]]
    })
    retire <- sum(vapply(offers, "[[", 0, "cost")) < z$price * (sum(p[k, ] > 0) - 1) - 1e-05
    for (i in seq_along(offers)) {
      f <- offers[[i]]
      if (!retire && f$cost >= -1e-05)
        next
      s <- members[i]
      x <- z$data[[s]]
      rec <- f$rec
      cc_set_h(x, f$r, c(k, rec), cbind(0, f$v))
      cc_set_mu(x, f$r, f$ion, f$m)
      m[s, k] <- 0
      m[s, rec] <- colSums(f$v)
      n <- n + 1L
      z$bands[rec, 1] <- pmin(z$bands[rec, 1], z$bands[k, 1])
      z$bands[rec, 2] <- pmax(z$bands[rec, 2], z$bands[k, 2])
    }
    if (!any(m[, k] > 0))
      cc_write_block(p, k, seq_len(ncol(p)), matrix(0, 1L, ncol(p)))
  }
  n
}
cc_mixture_trial <- function(bs, p, w, price, retain) {
  q <- p
  r <- if (retain)
    pmax(p[1, ] - drop(w %*% p[-1, , drop = FALSE]), 0) else numeric(ncol(p))
  size <- sum(r)
  q[1, ] <- r/max(size, 1e-80)
  vs <- lapply(bs, function(b) {
    v <- b$h
    for (j in seq_along(w)) if (w[j] > 0 && any(b$h[, 1] > 0))
      v[, j + 1L] <- cc_area_shape(b$t, b$h[, j + 1L] + w[j] * b$h[, 1])
    v[, 1] <- v[, 1] * size
    v
  })
  ms <- Map(function(b, v) b$off + v %*% q, bs, vs)
  best <- NULL
  for (step in 1:2) {
    for (i in seq_along(bs)) {
      f <- cc_coordinate_profiles(bs[[i]], vs[[i]], q, ms[[i]], price, only_present = TRUE)
      vs[[i]] <- f$v
      ms[[i]] <- f$m
    }
    fit <- cc_local_spectra(bs, vs, q, ms, price)
    vs <- fit$v
    q <- fit$p
    ms <- fit$m
    delta <- cc_trial_cost(bs, vs, q, ms, p, price)
    if (delta$cost < -1e-05 && (is.null(best) || delta$cost < best$cost))
      best <- c(fit, delta)
  }
  best
}
cc_shared_exchange <- function(z, full = FALSE) {
  p <- z$p
  m <- cc_mass(z$data)
  ns <- colSums(m > 0)
  b <- z$bands
  unit <- cc_unit(p)
  change <- 0
  if (full) {
    b <- cbind(rep(Inf, nrow(p)), rep(-Inf, nrow(p)))
    for (x in z$data) for (k in which(colSums(x$h) > 0)) {
      ix <- which(x$h[, k] > 0)
      b[k, ] <- c(min(b[k, 1], x$warp[ix[1]]), max(b[k, 2], x$warp[utils::tail(ix, 1L)]))
    }
  }
  representation <- function(k, ids = NULL) {
    seen <- m[, k] > 0
    if (!any(seen) || all(seen) || !any(p[k, ] > 0))
      return(list())
    if (is.null(ids)) {
      c <- which(b[, 1] <= b[k, 2] & b[, 2] >= b[k, 1] & ns > ns[k] & colSums(m[!seen, , drop = FALSE] >
        0) > 0)
      if (!length(c))
        return(list())
      ids <- if (full)
        list(c) else cc_choices(k, c, unit)
    }
    out <- list()
    for (c in ids) {
      if (!length(c) || any(rowSums(p[c, , drop = FALSE]) <= 0))
        next
      w <- cc_cone(t(p[c, , drop = FALSE]), p[k, ])
      pos <- w > 0
      if (!any(pos) || (!full && !all(pos)))
        next
      rel <- sqrt(sum((p[k, ] - drop(w %*% p[c, , drop = FALSE]))^2)/sum(p[k, ]^2))
      out[[length(out) + 1L]] <- list(k = k, rec = c[pos], w = w[pos], rel = rel)
    }
    out
  }
  offers <- unlist(lapply(seq_len(nrow(p)), representation), recursive = FALSE)
  if (!length(offers))
    return(0)
  used <- rep(FALSE, nrow(p))
  for (g in offers[order(vapply(offers, "[[", 0, "rel"))]) {
    k <- g$k
    if ((!full && used[k]) || g$rel > 0.32)
      next
    o <- representation(k, if (full)
      NULL else list(g$rec))
    if (!length(o))
      next
    f <- o[[1L]]
    rec <- f$rec
    w <- f$w
    rel <- f$rel
    if (full && rel <= 0.32)
      while (length(w) > 1L) {
        plausible <- rec[-which.min(w * sqrt(rowSums(p[rec, , drop = FALSE]^2)))]
        alt <- representation(k, list(plausible))
        if (!length(alt))
          break
        v <- alt[[1L]]
        total <- sum(m[, k]/vapply(z$data, function(x) x$phi, 0))
        if ((v$rel^2 - rel^2) * max(1, total) > z$price)
          break
        rec <- v$rec
        w <- v$w
        rel <- v$rel
      }
    if (!length(rec) || rel > 0.32)
      next
    ids <- c(k, rec)
    bs <- cc_blocks(z, ids)
    if (!length(bs))
      next
    fits <- Filter(Negate(is.null), lapply(if (full)
      c(FALSE, TRUE) else c(TRUE, FALSE), function(keep) cc_mixture_trial(bs, p[ids, , drop = FALSE], w, z$price,
      keep)))
    if (!length(fits))
      next
    fit <- fits[[which.min(vapply(fits, "[[", 0, "cost"))]]
    cc_commit_block(z, bs, ids, fit)
    change <- change + fit$cost
    used[k] <- TRUE
    if (full) {
      for (i in seq_along(bs)) m[bs[[i]]$si, ids] <- colSums(fit$v[[i]])
      ns <- colSums(m > 0)
    }
  }
  change
}
cc_overlap_groups <- function(h) {
  live <- which(colSums(h) > 0)
  if (length(live) < 2L)
    return(list())
  limits <- vapply(live, function(j) range(which(h[, j] >= 0.5 * max(h[, j]))), integer(2))
  o <- order(limits[1, ], limits[2, ], live)
  limits <- limits[, o, drop = FALSE]
  groups <- cumsum(c(TRUE, limits[1, -1] > utils::head(cummax(limits[2, ]), -1)))
  Filter(function(x) length(x) > 1L, split(live[o], groups))
}
cc_refit_overlaps <- function(z) {
  ns <- colSums(cc_mass(z$data) > 0)
  change <- 0
  for (x in z$data) for (ids in cc_overlap_groups(x$h)) {
    if (length(ids) < 2L)
      next
    ion <- which(colSums(z$p[ids, , drop = FALSE]) > 0)
    q <- z$p[ids, ion, drop = FALSE]
    all <- seq_along(ids)
    candidates <- Filter(function(j) {
      A <- t(q[-j, , drop = FALSE])
      w <- cc_cone(A, q[j, ])
      sum((A %*% w - q[j, ])^2) <= 0.32^2 * sum(q[j, ]^2)
    }, all)
    if (!length(candidates))
      next
    pos <- which(rowSums(x$h[, ids, drop = FALSE]) > 0)
    r <- seq.int(pos[1], utils::tail(pos, 1))
    v0 <- x$h[r, ids, drop = FALSE]
    old <- x$mu[r, ion, drop = FALSE]
    b <- list(t = x$t[r], y = x$y[r, ion, drop = FALSE], off = pmax(old - v0 %*% q, x$bg[r, ion,
      drop = FALSE]), phi = x$phi)
    oldface <- sum(cc_faces(b$t, v0))
    last <- (ns[ids] == 1L) * (rowSums(q > 0) - 1L)
    refit <- function(keep, start) {
      ps <- q[keep, , drop = FALSE]
      v <- start[, keep, drop = FALSE]
      v <- cc_scan_unmix(b$y, b$off, ps, v, matrix(TRUE, nrow(v), ncol(v)), x$phi, TRUE)
      for (j in seq_along(keep)) v[, j] <- cc_area_shape(b$t, v[, j])
      m <- b$off + v %*% ps
      best <- NULL
      for (it in 1:3) {
        a <- cc_coordinate_profiles(b, v, ps, m, z$price)
        v <- a$v
        m <- a$m
        full <- v0 * 0
        full[, keep] <- v
        val <- cc_diff(b$y, m, old, x$phi) + z$price * (sum(cc_faces(b$t, v)) - oldface - sum(last[colSums(full) ==
          0]))
        if (is.null(best) || val < best$cost)
          best <- list(cost = val, v = full, m = m, keep = keep)
      }
      best
    }
    active <- all
    best <- list(cost = 0, v = v0, m = old)
    for (depth in seq.int(0L, length(candidates))) {
      subsets <- if (depth == 0L)
        list(active) else list()
      if (length(active) > 1L)
        subsets <- c(subsets, lapply(intersect(active, candidates), function(j) setdiff(active,
          j)))
      if (!length(subsets))
        break
      opts <- lapply(subsets, refit, start = best$v)
      f <- opts[[which.min(vapply(opts, "[[", 0, "cost"))]]
      if (f$cost >= best$cost - 1e-05)
        break
      best <- f
      if (length(f$keep) == length(active))
        break
      active <- f$keep
    }
    if (best$cost < -1e-05) {
      cc_set_h(x, r, ids, best$v)
      cc_set_mu(x, r, ion, best$m)
      change <- change + best$cost
      retired <- ids[colSums(best$v) == 0]
      ns[retired] <- ns[retired] - 1L
    }
  }
  change
}
