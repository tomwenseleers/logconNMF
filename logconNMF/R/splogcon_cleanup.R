# Complete model moves on the original channels; likelihood kernels are native.
slcn_cleanup_delta <- function(f, Y, before, after, rows) {
  slcn_delta(Y[rows, , drop = FALSE], before, after, f$spec$kind,
    f$spec$dispersion, f$spec$theta, as.numeric(f$weights), nrow(Y),
    min(rows) - 1L, f$control$censor)
}

slcn_cleanup_prune <- function(f, Y) {
  audit <- list()
  calculate <- function(k) {
    rows <- which(f$H[, k] > 0)
    after <- pmax(f$mu[rows, , drop = FALSE] -
      as.numeric(f$H[rows, k]) %o% f$S[k, ], 0)
    delta <- slcn_cleanup_delta(f, Y, f$mu[rows, , drop = FALSE], after, rows) - f$lambda
    if (is.finite(delta)) delta else Inf
  }
  deltas <- vapply(seq_len(ncol(f$H)), calculate, 0)
  while (length(deltas) && min(deltas) < -1e-6) {
    k <- which.min(deltas);rows <- which(f$H[, k] > 0)
    audit[[length(audit) + 1L]] <- data.frame(component = paste(f$members[[k]], collapse = ";"), delta = deltas[k])
    f$mu[rows, ] <- pmax(f$mu[rows, , drop = FALSE] - as.numeric(f$H[rows, k]) %o% f$S[k, ], 0)
    f$H <- f$H[, -k, drop = FALSE];f$S <- f$S[-k, , drop = FALSE]
    f$members <- f$members[-k];deltas <- deltas[-k]
    # Disjoint temporal supports have unchanged conditional deletion costs.
    affected <- which(colSums(f$H[rows, , drop = FALSE]) > 0)
    deltas[affected] <- vapply(affected, calculate, 0)
  }
  list(fit = f, audit = if (length(audit)) do.call(rbind, audit) else data.frame())
}

slcn_cleanup_merge <- function(f, Y, group) {
  keep <- setdiff(seq_len(ncol(f$H)), group)
  raw <- as.numeric(f$H[, group, drop = FALSE] %*% rowSums(f$S[group, , drop = FALSE]))
  at <- which(raw > max(raw) * 1e-10)
  rows <- seq.int(max(1L, min(at) - 3L), min(nrow(Y), max(at) + 3L))
  fixed <- pmax(f$mu[rows, , drop = FALSE] -
    as.matrix(f$H[rows, group, drop = FALSE] %*% f$S[group, , drop = FALSE]), 0)
  h <- spolo_lc(raw[rows])$h
  h <- pmax(h, max(h) * 1e-12);h <- h / sum(h)
  s <- matrix(colSums(f$S[group, , drop = FALSE]), 1L)
  for (cycle in 1:2) {
    z <- slcn_joint_spectra(f, Y, matrix(h, ncol = 1L), s, fixed, rows)
    s <- z$S
    z <- slcn_log_profile(Y[rows, , drop = FALSE], z$mu, h, as.numeric(s),
      f$spec$kind, f$spec$dispersion, f$spec$theta, as.numeric(f$weights),
      nrow(Y), min(rows) - 1L, 4L, f$control$censor, 12L, as.integer(f$ncores), mode = 0L)
    h <- z$profile;scale <- sum(h)
    if (!is.finite(scale) || scale <= 0) return(NULL)
    h <- h / scale;s <- s * scale
  }
  z <- slcn_joint_spectra(f, Y, matrix(h, ncol = 1L), s, fixed, rows)
  delta <- slcn_cleanup_delta(f, Y, f$mu[rows, , drop = FALSE], z$mu, rows) -
    f$lambda * (length(group) - 1L)
  if (!is.finite(delta) || delta >= -1e-6) return(NULL)
  f$H <- cbind(f$H[, keep, drop = FALSE], sparseMatrix(i = rows,
    j = rep.int(1L, length(rows)), x = h, dims = c(nrow(Y), 1L)))
  f$S <- rbind(f$S[keep, , drop = FALSE], z$S);f$mu[rows, ] <- z$mu
  f$members <- c(f$members[keep], list(sort(unlist(f$members[group]))))
  list(fit = f, delta = delta)
}

slcn_cleanup_fuse <- function(f, Y, max_merges = 12L, tries = 20L) {
  audit <- list()
  for (pass in seq_len(max_merges)) {
    if (ncol(f$H) < 2L) break
    ch <- as.matrix(crossprod(f$H)) / sqrt(colSums(f$H^2) %o% colSums(f$H^2))
    cs <- tcrossprod(f$S) / sqrt(rowSums(f$S^2) %o% rowSums(f$S^2))
    pairs <- which(upper.tri(ch) & is.finite(cs) &
      (ch >= .95 | (cs >= .995 & ch > .1)), arr.ind = TRUE)
    if (!nrow(pairs)) break
    order <- order(ch[pairs] + pmax(cs[pairs] - .95, 0), decreasing = TRUE)
    best <- NULL;delta <- -1e-6
    for (i in head(order, tries)) {
      z <- slcn_cleanup_merge(f, Y, as.integer(pairs[i, ]))
      if (!is.null(z) && z$delta < delta) {
        best <- z;delta <- z$delta;pair <- pairs[i, ]
      }
    }
    if (is.null(best)) break
    audit[[length(audit) + 1L]] <- data.frame(
      a = paste(f$members[[pair[1]]], collapse = ";"),
      b = paste(f$members[[pair[2]]], collapse = ";"), delta = delta)
    f <- best$fit
  }
  list(fit = f, audit = if (length(audit)) do.call(rbind, audit) else data.frame())
}

slcn_cleanup_residual <- function(f, Y, X, width = 512L, candidates = 16L, max_add = 6L) {
  audit <- list();n <- nrow(Y);m <- ncol(Y);tested <- integer()
  for (left in seq.int(1L, n, by = width)) {
    window <- seq.int(left, min(n, left + width - 1L))
    columns <- which(colSums(X[window, , drop = FALSE]) > .5)
    if (!length(columns)) next
    rows <- which(rowSums(X[, columns, drop = FALSE]) > 0)
    rows <- seq.int(min(rows), max(rows));D <- X[rows, columns, drop = FALSE]
    mu <- pmax(f$mu[rows, , drop = FALSE], .1)
    variance <- switch(f$spec$name, gaussian = matrix(f$spec$dispersion, nrow(mu), m),
      poisson = mu, quasipoisson = mu * f$spec$dispersion,
      nb = mu + mu^2 / f$spec$theta, nb1 = mu * f$spec$dispersion,
      Gamma_scale = mu * f$spec$dispersion, Gamma = mu^2 * f$spec$dispersion)
    W <- if (length(f$weights) == 1L) f$weights else f$weights[rows, , drop = FALSE]
    G <- as.matrix(crossprod(D, W * (Y[rows, , drop = FALSE] - mu) / variance))
    I <- as.matrix(crossprod(D^2, W / variance))
    gain <- rowSums(pmax(G, 0)^2 / pmax(I, 1e-20))
    # A working-score screen proposes moves; the actual family likelihood decides them.
    added <- attempted <- 0L
    for (j in order(gain, decreasing = TRUE)) {
      pos <- columns[j]
      if (any(abs(pos - tested) < 3L)) next
      tested <- c(tested, pos);attempted <- attempted + 1L
      at <- which(X[, pos] > 0);h <- as.numeric(X[at, pos]);h <- h / sum(h)
      fixed <- f$mu[at, , drop = FALSE]
      s <- matrix(pmax(as.numeric(crossprod(h, pmax(Y[at, , drop = FALSE] - fixed, 0))), 1e-6), 1L)
      z <- slcn_joint_spectra(f, Y, matrix(h, ncol = 1L), s, fixed, at)
      delta <- slcn_cleanup_delta(f, Y, fixed, z$mu, at) + f$lambda
      ok <- is.finite(delta) && delta < -1e-6
      audit[[length(audit) + 1L]] <- data.frame(scan = pos, score = gain[j], delta = delta, accepted = ok)
      if (ok) {
        f$H <- cbind(f$H, sparseMatrix(i = at, j = rep.int(1L, length(at)), x = h, dims = c(n, 1L)))
        f$S <- rbind(f$S, z$S);f$mu[at, ] <- z$mu
        f$members <- c(f$members, list(-as.integer(pos)))
        proposal <- slcn_refine_block(f, Y, block = ncol(f$H), cycles = 2L,
          profile_steps = 3L, adaptive = FALSE, joint = TRUE)
        # Compare every row touched by the local profile/background routine.
        dd <- slcn_cleanup_delta(f, Y, f$mu, proposal$mu, seq_len(n)) +
          f$lambda * (ncol(proposal$H) - ncol(f$H))
        if (is.finite(dd) && dd <= 0) f <- proposal
        added <- added + 1L
      }
      if (attempted >= candidates || added >= max_add) break
    }
  }
  list(fit = f, audit = if (length(audit)) do.call(rbind, audit) else data.frame())
}

slcn_cleanup <- function(f, Y, X, verbose = FALSE) {
  Y <- as.matrix(Y);storage.mode(Y) <- "double"
  if (is.null(f$mu)) f$mu <- as.matrix(f$H %*% f$S) + f$background %*% f$alpha + f$offset
  if (length(f$members) != ncol(f$H)) f$members <- lapply(seq_len(ncol(f$H)), identity)
  criterion <- function(z) -2 * slcn_loglik(Y, z$mu, z$spec$kind, z$spec$dispersion,
    z$spec$theta, as.numeric(z$weights), as.integer(z$ncores), z$control$censor) + z$lambda * ncol(z$H)
  initial <- criterion(f);audit <- list()
  stage <- function(name, expression) {
    before <- criterion(f);seconds <- slcn_time(z <- force(expression))[["elapsed"]]
    delta <- criterion(z$fit) - before;ok <- is.finite(delta) && delta <= 1e-6
    if (ok) f <<- z$fit
    audit[[name]] <<- z$audit
    f$history <<- rbind(f$history, data.frame(stage = name, seconds = seconds,
      delta = delta, accepted = ok, compounds = ncol(z$fit$H)))
    if (verbose) {cat(name, round(seconds, 2), "s; change", signif(delta, 7),
      if (ok) "accepted;" else "rejected;", ncol(f$H), "compounds\n");flush.console()}
  }
  stage("original_prune", slcn_cleanup_prune(f, Y))
  stage("original_fusion", slcn_cleanup_fuse(f, Y))
  stage("original_residual", slcn_cleanup_residual(f, Y, X))
  # Reconstruct once to remove roundoff accumulated by local subtraction/addition.
  f$mu <- as.matrix(f$H %*% f$S) + f$background %*% f$alpha + f$offset
  f$loglik <- slcn_loglik(Y, f$mu, f$spec$kind, f$spec$dispersion, f$spec$theta,
    as.numeric(f$weights), as.integer(f$ncores), f$control$censor)
  f$rank <- ncol(f$H);f$ic_df <- ncol(Y) * f$rank + f$nuisance_df
  f$objective <- -2 * f$loglik + f$lambda * f$rank
  f$complete_delta <- f$objective - f$initial_objective
  f$ledger_error <- f$complete_delta - sum(f$history$delta[f$history$accepted])
  f$cleanup <- list(initial_objective = initial, delta = f$objective - initial, moves = audit)
  f$state_ids <- f$source_ids <- NULL
  slcn_temporal_order(f)
}
