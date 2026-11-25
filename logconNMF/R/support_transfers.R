# Screen removals at four times the fitted penalty. Only the ordinary
# original-ion objective can accept a change in a fitted occurrence or identity.
cc_profile_ranges <- function(h) {
  a <- as(h, "dgCMatrix")
  live <- which(diff(a@p) > 0L)
  bounds <- cbind(rep(Inf, ncol(h)), rep(-Inf, ncol(h)))
  bounds[live, 1L] <- a@i[1L + a@p[live]] + 1L
  bounds[live, 2L] <- a@i[a@p[live + 1L]] + 1L
  bounds
}

cc_screen_support <- function(z) {
  p <- z$p
  rank <- nrow(p)
  samples <- length(z$data)
  first <- second <- matrix(0L, rank, samples)
  fraction <- matrix(0, rank, samples)
  ranges <- lapply(z$data, function(x) cc_profile_ranges(x$h))
  norm2 <- rowSums(p * p)
  # Bound the dense Gram allocation. Large dictionaries need inner products
  # only for spectra whose fitted profiles intersect in the current sample.
  gram <- if (as.double(rank) * rank <= 2^22) tcrossprod(p)
  dot <- function(k, j) {
    if (is.null(gram)) return(drop(p[j, , drop = FALSE] %*% p[k, ]))
    gram[j, k]
  }
  for (s in seq_along(ranges)) {
    b <- ranges[[s]]
    for (k in which(is.finite(b[, 1L]))) {
      j <- setdiff(which(b[, 1L] <= b[k, 2L] & b[, 2L] >= b[k, 1L]), k)
      if (!length(j)) next
      inner <- dot(k, j)
      cosine <- inner / sqrt(pmax(norm2[k] * norm2[j], 1e-80))
      closest <- which.max(cosine)
      if (cosine[closest] <= 0) next
      a <- first[k, s] <- j[closest]
      other <- which(cosine > 0 & j != a)
      if (!length(other)) next
      candidates <- j[other]
      cross <- dot(a, candidates)
      length2 <- norm2[candidates] + norm2[a] - 2 * cross
      projection <- inner[other] - inner[closest] - cross + norm2[a]
      weight <- pmin(1, pmax(0, projection / pmax(length2, 1e-80)))
      gain <- 2 * weight * projection - weight^2 * length2
      at <- which.max(gain)
      if (gain[at] > 0) {
        second[k, s] <- candidates[at]
        fraction[k, s] <- weight[at]
      }
    }
  }
  cost <- cc_screen_transfers(z$data, p, first, second, fraction, z$threads)
  faces <- matrix(unlist(lapply(z$data, function(x) cc_faces(x$t, x$h)),
                         use.names = FALSE), rank, samples)
  spectral <- pmax(0, rowSums(p > 0) - 1L)
  excess <- cost - 4 * z$price * faces
  priority <- rowSums(pmax(excess, 0), na.rm = TRUE) - 4 * z$price * spectral
  shared <- priority <= 0 & rowSums(!is.na(cost)) > 0
  list(flagged = !is.na(cost) & (excess < 0 | shared), shared = shared,
       priority = priority, spectral = spectral, first = first, second = second,
       fraction = fraction, ranges = ranges)
}

# Remove the first profile and allow at most two overlapping profiles to
# absorb its signal. Their spectra stay fixed; everything else is in the offset.
cc_transfer_proposal <- function(x, p, ids, bounds, price, weights = 1) {
  rows <- seq.int(min(bounds[ids, 1L]), max(bounds[ids, 2L]))
  ions <- which(colSums(p[ids, , drop = FALSE]) > 0)
  q <- p[ids, ions, drop = FALSE]
  old <- x$mu[rows, ions, drop = FALSE]
  h <- x$h[rows, ids, drop = FALSE]
  block <- list(y = x$y[rows, ions, drop = FALSE], t = x$t[rows], phi = x$phi,
                off = pmax(old - h %*% q, x$bg[rows, ions, drop = FALSE]))
  old_df <- sum(cc_faces(block$t, h))
  score <- function(v, mean) cc_diff(block$y, mean, old, block$phi) +
    price * (sum(cc_faces(block$t, v)) - old_df)
  v <- h
  v[, 1L] <- 0
  mean <- block$off + v %*% q
  best <- list(rows = rows, ions = ions, ids = ids, h = v, mu = mean, delta = score(v, mean))
  partners <- seq_along(ids)[-1L]
  if (!length(partners)) return(best)
  profiles <- h[, partners, drop = FALSE]
  for (j in seq_along(partners))
    profiles[, j] <- cc_area_shape(block$t, profiles[, j] + weights[j] * h[, 1L])
  spectra <- q[partners, , drop = FALSE]
  mean <- block$off + profiles %*% spectra
  for (i in seq_len(2L)) {
    update <- cc_coordinate_profiles(block, profiles, spectra, mean, price)
    profiles <- update$v
    mean <- update$m
    if (any(colSums(profiles) <= 0)) break
    v[, partners] <- profiles
    delta <- score(v, mean)
    if (delta < best$delta) {
      best$h <- v
      best$mu <- mean
      best$delta <- delta
    } else if (length(partners) == 1L) break
  }
  best
}

cc_try_transfers <- function(z, scan) {
  ranges <- scan$ranges
  unit <- cc_unit(z$p)
  candidates <- which(rowSums(scan$flagged) > 0)
  removed <- 0L
  for (k in candidates[order(scan$priority[candidates])]) {
    live <- which(vapply(ranges, function(b) is.finite(b[k, 1L]), TRUE))
    samples <- live[scan$flagged[k, live]]
    if (!length(samples)) next
    proposals <- lapply(samples, function(s) {
      b <- ranges[[s]]
      neighbours <- setdiff(which(b[, 1L] <= b[k, 2L] & b[, 2L] >= b[k, 1L]), k)
      closest <- integer()
      if (length(neighbours)) {
        cosine <- drop(unit[neighbours, , drop = FALSE] %*% unit[k, ])
        if (max(cosine) > 0) closest <- neighbours[which.max(cosine)]
      }
      trial <- cc_transfer_proposal(z$data[[s]], z$p, c(k, closest), b, z$price)
      pair <- c(scan$first[k, s], scan$second[k, s])
      if (all(pair %in% neighbours)) {
        w <- scan$fraction[k, s]
        both <- cc_transfer_proposal(z$data[[s]], z$p, c(k, pair), b, z$price, c(1 - w, w))
        if (both$delta < trial$delta) trial <- both
      }
      trial
    })
    delta <- vapply(proposals, function(u) u$delta, 0)
    spectral_saving <- z$price * scan$spectral[k]
    retire <- length(live) == length(samples) && sum(pmax(delta, 0)) <= spectral_saving
    take <- if (retire) seq_along(samples) else which(delta < -1e-6)
    total <- sum(delta[take]) - if (retire) spectral_saving else 0
    if (!length(take) || total >= -1e-6) next
    for (at in take) {
      s <- samples[at]
      u <- proposals[[at]]
      cc_set_mu(z$data[[s]], u$rows, u$ions, u$mu)
      cc_set_h(z$data[[s]], u$rows, u$ids, u$h)
      ranges[[s]][k, ] <- c(Inf, -Inf)
      for (j in seq_along(u$ids)[-1L])
        ranges[[s]][u$ids[j], ] <- range(u$rows[u$h[, j] > 0])
    }
    removed <- removed + length(take)
  }
  removed
}
