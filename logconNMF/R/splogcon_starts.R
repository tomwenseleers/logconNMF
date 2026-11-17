# Aggregated proposals are used only under a closed Poisson working model.
# Other families, observation weights and offsets retain full-channel screening.
slcn_restart_groups <- function(f, Y) {
  group <- f$aggregation$group
  if (is.null(group) || length(group) != ncol(Y) ||
      !f$spec$kind %in% 1:2 || length(f$spec$dispersion) != 1L ||
      length(f$weights) != 1L || any(f$offset != 0)) return(NULL)
  group
}

# Alternative local starts for nearly coincident component modes. Screening
# scores never authorize a move: completed candidates and controls use Y.
slcn_profile_starts <- function(f, Y, cycles = NULL, knots = 128L) {
  begun <- proc.time()[["elapsed"]]
  f <- slcn_sync(f, Y);K <- ncol(f$H)
  if (K < 2L) return(f)
  geometry <- slcn_profile_geometry(f)
  sorted <- order(geometry$apex);eligible <- rep(FALSE, K)
  for (a in seq_len(K - 1L)) {
    k <- sorted[a];b <- a + 1L
    while (b <= K && geometry$apex[sorted[b]] - geometry$apex[k] <= geometry$width[k]/4) {
      j <- sorted[b]
      if (geometry$apex[j] - geometry$apex[k] <= geometry$width[j]/4) eligible[c(k,j)] <- TRUE
      b <- b + 1L
    }
  }
  if (!any(eligible)) return(f)
  group <- slcn_restart_groups(f, Y)
  aggregated <- !is.null(group)
  if (is.null(cycles)) cycles <- if (aggregated) 8L else 16L
  screen_Y <- Y
  compress <- identity
  expand <- function(z, original, coarse) z
  if (aggregated) {
    P <- sparseMatrix(i = seq_along(group), j = group, x = 1,
      dims = c(ncol(Y), max(group)))
    screen_Y <- as.matrix(Y %*% P)
    compress <- function(z) {
      z$S <- as.matrix(z$S %*% P);z$alpha <- as.matrix(z$alpha %*% P)
      z$mu <- as.matrix(z$mu %*% P);z$y <- screen_Y
      z$offset <- 0;z$aggregation <- NULL
      z
    }
    expand <- function(z, original, coarse) {
      fractions <- original$S / pmax(coarse$S[, group, drop = FALSE], 1e-30)
      original$H <- z$H
      original$S <- z$S[, group, drop = FALSE] * fractions
      # Restart blocks keep nuisance terms fixed; do not round-trip alpha.
      original$mu <- as.matrix(original$H %*% original$S) +
        original$background %*% original$alpha + original$offset
      original
    }
  }
  criterion <- function(z, data) slcn_loglik(data, z$mu, z$spec$kind,
    z$spec$dispersion, z$spec$theta, as.numeric(z$weights),
    as.integer(z$ncores), z$control$censor)
  deform <- function(h, left, right) {
    at <- which(h > 0);apex <- which.max(h);grid <- seq_along(h)
    position <- apex + (grid - apex) / ifelse(grid < apex, left, right)
    h <- exp(stats::approx(at, log(h[at]), xout = position,
      yleft = -Inf, yright = -Inf)$y)
    h / sum(h)
  }
  coarse <- compress(f);cache <- new.env(parent = emptyenv())
  hits <- misses <- 0L;audit <- list();sources <- f$members
  for (k in which(eligible)[order(rowSums(f$S)[eligible], decreasing = TRUE)]) {
    h <- as.numeric(f$H[, k])
    if (sum(h > 0) < 3L) next
    geometry <- slcn_profile_geometry(f)
    trials <- lapply(list(c(2,1), c(1,2), c(2,2)),
      function(scale) deform(h, scale[1L], scale[2L]))
    at <- which(Reduce(`+`, trials) + h > max(h) * 1e-6)
    block <- which(geometry$left <= max(at) & geometry$right >= min(at))
    screen <- function(z) slcn_refine_block(z, screen_Y, block = block,
      cycles = 1L, profile_steps = 2L, adaptive = FALSE, joint = TRUE,
      pad = 2L, free = block, knots = 16L)
    key <- paste(block, collapse = ",")
    if (exists(key, cache, inherits = FALSE)) {
      control_score <- get(key, cache);hits <- hits + 1L
    } else {
      control <- screen(coarse)
      control_score <- if (ncol(control$H) == K) criterion(control, screen_Y) else -Inf
      # Only this scalar is reused; caching models would retain full mean arrays.
      assign(key, control_score, cache);misses <- misses + 1L
    }
    before <- criterion(coarse, screen_Y)
    best <- NULL;best_score <- before;choice <- "none"
    for (j in seq_along(trials)) {
      z <- coarse;z$H[, k] <- trials[[j]]
      z$mu <- as.matrix(z$H %*% z$S) + z$background %*% z$alpha + z$offset
      z <- screen(z)
      score <- if (ncol(z$H) == K) criterion(z, screen_Y) else -Inf
      if (is.finite(score) && score > max(best_score, control_score) + 1e-6) {
        best <- z;best_score <- score;choice <- c("left", "right", "both")[j]
      }
    }
    accepted <- FALSE;delta <- 0
    if (!is.null(best)) {
      current <- criterion(f, Y)
      polish <- function(z, iterations = cycles) slcn_refine_block(z, Y,
        block = block, cycles = iterations, profile_steps = 6L,
        adaptive = FALSE, joint = TRUE, pad = 2L, free = block, knots = knots)
      candidate <- polish(expand(best, f, coarse))
      control <- polish(f, cycles + 1L)
      score <- if (ncol(candidate$H) == K) criterion(candidate, Y) else -Inf
      reference <- if (ncol(control$H) == K) criterion(control, Y) else -Inf
      if (is.finite(score) && score > max(current, reference) + 1e-6) {
        delta <- -2 * (score - current);f <- candidate;coarse <- compress(f)
        accepted <- TRUE;rm(list = ls(cache), envir = cache)
      }
    }
    audit[[length(audit) + 1L]] <- data.frame(component = k, start = choice,
      screen_delta = -2 * (best_score - before), delta = delta, accepted = accepted)
  }
  f <- slcn_sync(f, Y)
  f$profile_starts <- list(audit = if (length(audit)) do.call(rbind, audit) else data.frame(),
    seconds = proc.time()[["elapsed"]] - begun, source_components = sources,
    screen_channels = ncol(screen_Y), aggregated = aggregated,
    polish_cycles = cycles, cache_hits = hits, cache_misses = misses)
  f
}
