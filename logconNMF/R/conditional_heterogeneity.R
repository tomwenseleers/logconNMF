cc_partition <- function(bs, p) {
  if (length(bs) < 2L)
    return(NULL)
  W <- lapply(bs, function(b) 1/(b$old * b$phi))
  g <- do.call(rbind, Map(function(b, w) colSums((b$y - b$old) * w * as.numeric(b$h)), bs, W))
  I <- pmax(do.call(rbind, Map(function(b, w) colSums(w * (as.numeric(b$h)^2)), bs, W)), 1e-80)
  s <- drop(g %*% p)/pmax(drop(I %*% (p * p)), 1e-80)
  g <- g - s * sweep(I, 2L, p, "*")
  f <- pmax(sweep(g/I, 2L, p, "+"), 0)
  r <- sqrt(f/pmax(rowSums(f), 1e-60))
  v <- svd(sweep(r, 2L, colMeans(r)), nu = 0L, nv = 1L)$v[, 1L]
  order <- order(drop(r %*% v))
  n <- length(bs)
  cuts <- c(lapply(seq_len(n - 1L), function(i) order[seq_len(i)]), as.list(seq_len(n)))
  gain <- function(a, b) {
    d <- pmax(-p, a/b)
    sum(2 * a * d - b * d * d)
  }
  base <- gain(colSums(g), colSums(I))
  best <- NULL
  for (inc in cuts) {
    pos <- seq_len(n) %in% inc
    G <- gain(colSums(g[pos, , drop = FALSE]), colSums(I[pos, , drop = FALSE])) + gain(colSums(g[!pos,
      , drop = FALSE]), colSums(I[!pos, , drop = FALSE])) - base
    if (is.null(best) || G > best$gain)
      best <- list(gain = G, mask = pos)
  }
  best
}
cc_train_ray <- function(bs, p, price, j = 1L) {
  data <- function(key) do.call(rbind, lapply(bs, function(b) b[[key]]/b$phi))
  cc_local_spectrum(data("y"), data("off"), unlist(lapply(bs, function(b) b$h[, j]/b$phi), use.names = FALSE),
    p, price)
}
cc_replicated <- function(bs, p, mask, price) {
  if (min(sum(mask), sum(!mask)) < 2L)
    return(FALSE)
  for (grp in list(which(mask), which(!mask))) {
    gain <- 0
    for (i in grp) {
      b <- bs[[i]]
      local <- cc_train_ray(bs[setdiff(grp, i)], p, price)
      pooled <- cc_train_ray(bs[-i], p, price)
      if (sum(local) == 0 || sum(pooled) == 0)
        return(FALSE)
      local <- local/sum(local)
      pooled <- pooled/sum(pooled)
      calc <- function(q) {
        h <- cc_local_profile(b$t, b$y, b$off, q, b$h[, 1], b$phi, price, TRUE)
        b$off + tcrossprod(h, q)
      }
      gain <- gain + cc_diff(b$y, calc(pooled), calc(local), b$phi)
    }
    if (gain <= 0)
      return(FALSE)
  }
  TRUE
}
cc_split_offer <- function(bs, p, mask, price) {
  out <- list()
  delta <- 0
  df <- 1 - sum(p > 0)
  for (g in list(which(mask), which(!mask))) {
    q <- cc_train_ray(bs[g], p, price)
    s <- sum(q)
    if (s <= 0)
      return(NULL)
    q <- q/s
    ion <- which(q > 0)
    df <- df + length(ion) - 1L
    vals <- lapply(bs[g], function(b) {
      v <- cc_local_profile(b$t, b$y[, ion, drop = FALSE], b$off[, ion, drop = FALSE], q[ion],
        b$h[, 1] * s, b$phi, price)
      m <- b$off + tcrossprod(v, q)
      list(v = v, m = m)
    })
    if (!any(vapply(vals, function(v) sum(v$v) > 0, TRUE)))
      return(NULL)
    for (i in seq_along(g)) {
      b <- bs[[g[i]]]
      f <- vals[[i]]
      delta <- delta + cc_diff(b$y, f$m, b$old, b$phi)
      df <- df + cc_faces(b$t, as.matrix(f$v)) - cc_faces(b$t, b$h)
    }
    out[[length(out) + 1L]] <- list(p = q, values = vals, at = g)
  }
  if (delta + price * df >= -1e-05)
    NULL else list(groups = out, cost = delta + price * df)
}
cc_two_ray <- function(b, p, v) {
  for (it in 1:4) {
    w <- 1/pmax(b$off + v %*% p, 1e-08)
    a <- drop(w %*% (p[1, ]^2))
    cross <- drop(w %*% (p[1, ] * p[2, ]))
    c <- drop(w %*% (p[2, ]^2))
    u <- ((b$y - b$off) * w) %*% t(p)
    det <- a * c - cross^2
    stable <- det > 1e-12 * a * c
    v[, ] <- 0
    v[stable, 1] <- (c * u[, 1] - cross * u[, 2])[stable]/det[stable]
    v[stable, 2] <- (a * u[, 2] - cross * u[, 1])[stable]/det[stable]
    edge <- which(!stable | v[, 1] < 0 | v[, 2] < 0)
    if (length(edge)) {
      s <- pmax(u, 0)/pmax(cbind(a, c), 1e-100)
      sel <- max.col(s * u, ties.method = "first")
      v[edge, ] <- 0
      at <- cbind(edge, sel[edge])
      v[at] <- s[at]
    }
  }
  v
}
cc_polish_pair <- function(z, ids) {
  bs <- cc_blocks(z, ids)
  seen <- do.call(rbind, lapply(bs, function(b) colSums(b$h) > 0))
  if (!any(rowSums(seen) == 2L) || !all(colSums(seen[rowSums(seen) == 1L, , drop = FALSE]) > 0))
    return(0)
  p <- oldp <- z$p[ids, , drop = FALSE]
  for (j in 1:2) {
    q <- cc_train_ray(bs[seen[, j] & !seen[, 3L - j]], oldp[j, ], z$price, j)
    if (sum(q) <= 0)
      return(0)
    p[j, ] <- q/sum(q)
  }
  vs <- lapply(bs, "[[", "h")
  best <- NULL
  for (it in 1:2) {
    ms <- vector("list", length(bs))
    for (i in seq_along(bs)) {
      b <- bs[[i]]
      v <- vs[[i]]
      js <- which(seen[i, ])
      if (length(js) == 2L) {
        v <- cc_two_ray(b, p, v)
        for (j in js) v[, j] <- cc_area_shape(b$t, v[, j])
      }
      upd <- cc_coordinate_profiles(b, v, p, b$off + v %*% p, z$price, js)
      vs[[i]] <- upd$v
      ms[[i]] <- upd$m
    }
    fit <- cc_local_spectra(bs, vs, p, ms, z$price)
    vs <- fit$v
    p <- fit$p
    ms <- fit$m
    if (any(rowSums(p) == 0))
      break
    dh <- cc_trial_cost(bs, vs, p, ms, oldp, z$price)
    if (dh$cost < -1e-05 && (is.null(best) || dh$cost < best$cost))
      best <- c(fit, dh)
  }
  if (is.null(best))
    return(0)
  cc_commit_block(z, bs, ids, best)
  best$cost
}
cc_split_spectra <- function(z) {
  p <- z$p
  op <- numeric(nrow(p))
  change <- 0
  for (k in seq_len(nrow(p))) if (any(p[k, ] > 0)) {
    part <- cc_partition(cc_blocks(z, k), p[k, ])
    if (!is.null(part))
      op[k] <- part$gain - z$price * (sum(p[k, ] > 0) - 1)
  }
  for (k in order(-op)) {
    if (op[k] <= 0)
      break
    p <- z$p[k, ]
    bs <- cc_blocks(z, k)
    part <- cc_partition(bs, p)
    if (is.null(part) || part$gain <= z$price * (sum(p > 0) - 1) || !cc_replicated(bs, p, part$mask,
      z$price))
      next
    offer <- cc_split_offer(bs, p, part$mask, z$price)
    if (is.null(offer))
      next
    groups <- offer$groups
    change <- change + offer$cost
    oldrank <- nrow(z$p)
    j <- oldrank + 1L
    ids <- c(k, j)
    z$p <- rbind(z$p, groups[[2]]$p)
    z$p[k, ] <- groups[[1]]$p
    z$bands <- rbind(z$bands, z$bands[k, ])
    for (x in z$data) x$h <- cbind(x$h, numeric(length(x$t)))
    for (l in 1:2) {
      s <- groups[[l]]
      for (i in seq_along(s$at)) {
        b <- bs[[s$at[i]]]
        x <- z$data[[b$si]]
        v <- s$values[[i]]
        cc_set_h(x, b$r, k, numeric(length(b$r)))
        cc_set_h(x, b$r, ids[l], v$v)
        cc_set_mu(x, b$r, seq_along(p), v$m)
      }
    }
    ion <- which(colSums(z$p[ids, , drop = FALSE]) > 0)
    ps <- z$p[ids, ion, drop = FALSE]
    coexisting <- FALSE
    for (b in bs) {
      x <- z$data[[b$si]]
      prior <- x$h[b$r, ids, drop = FALSE]
      old <- x$mu[b$r, ion, drop = FALSE]
      b$y <- b$y[, ion, drop = FALSE]
      b$off <- b$off[, ion, drop = FALSE]
      v <- cc_scan_unmix(b$y, b$off, ps, prior, matrix(TRUE, nrow(prior), 2L), b$phi)
      for (l in 1:2) v[, l] <- cc_area_shape(b$t, v[, l])
      f <- cc_coordinate_profiles(b, v, ps, b$off + v %*% ps, z$price)
      gain <- cc_diff(b$y, f$m, old, b$phi) + z$price * sum(cc_faces(b$t, f$v) - cc_faces(b$t,
        prior))
      if (gain < -1e-05) {
        cc_set_mu(x, b$r, ion, f$m)
        cc_set_h(x, b$r, ids, f$v)
        coexisting <- TRUE
        change <- change + gain
      }
    }
    if (coexisting)
      change <- change + cc_polish_pair(z, ids)
  }
  change
}
