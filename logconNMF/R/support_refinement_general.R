# These solvers already have family-aware likelihood, exposure and weight
# handling. Keep their original dimension conventions throughout the search.
slcn_end_refine <- function(f, y, control) {
  tic <- proc.time()[["elapsed"]]
  input <- f
  sweep <- function(f, ids, i) {
    next_fit <- slcn_refine_block(f, y, block = ids, cycles = 1L, adaptive = FALSE,
                                  joint = TRUE, background = length(ids) == ncol(f$H))
    next_fit <- slcn_sync(next_fit, y)
    if (is.finite(next_fit$objective) && next_fit$objective <= f$objective) next_fit else f
  }
  changes <- function(old, new) {
    ds <- rowSums(abs(old$S - new$S)) / pmax(rowSums(old$S), rowSums(new$S), 1e-12)
    pmax(ds, nmf_profile_change(old$H, new$H))
  }
  run <- nmf_sweep_schedule(f, seq_len(ncol(f$H)), control$refinement_sweeps, control$selective,
                             sweep, identity, changes)
  f <- run$state
  retained <- if (control$refinement_sweeps) "ordinary sweeps" else "input"
  if (control$penalty_search && f$lambda > 0 && ncol(f$H)) {
    candidate <- f
    candidate$members <- as.list(seq_len(ncol(f$H)))
    candidate$lambda <- 4 * f$lambda
    candidate <- slcn_cleanup_prune(candidate, y)$fit
    if (ncol(candidate$H) < ncol(f$H)) {
      keep <- unlist(candidate$members, use.names = FALSE)
      candidate$source_ids <- f$source_ids[keep]
      candidate$state_ids <- f$state_ids[keep]
      candidate$members <- NULL
      candidate$lambda <- f$lambda
      candidate <- slcn_sync(candidate, y)
      candidate <- sweep(candidate, seq_len(ncol(candidate$H)), 1L)
      if (candidate$objective < f$objective) {
        f <- candidate
        retained <- "fourfold deletion"
      }
    }
  }
  elapsed <- proc.time()[["elapsed"]] - tic
  f$refinement <- list(before = input$objective, after = f$objective, retained = retained,
                        visited = run$visited, possible = run$possible)
  f$history <- rbind(f$history, data.frame(stage = "terminal refinement", seconds = elapsed,
    delta = f$objective - input$objective, accepted = TRUE, compounds = ncol(f$H)))
  f$complete_delta <- f$objective - f$initial_objective
  f$ledger_error <- f$complete_delta - sum(f$history$delta[f$history$accepted])
  f$seconds["refinement"] <- f$seconds["refinement"] + elapsed
  colnames(f$S) <- colnames(y)
  dimnames(f$mu) <- dimnames(y)
  slcn_temporal_order(f)
}

sinmf_factor_changes <- function(old, new, ids) {
  change <- rowSums(abs(new$P[ids,,drop=FALSE] - old$P[ids,,drop=FALSE]))
  for (s in seq_along(old$paths)) {
    a <- readRDS(old$paths[s])
    b <- readRDS(new$paths[s])
    root <- old$map[a$gid]
    at <- match(root, ids)
    amplitudes <- function(f, scale) f$H %*% Matrix::Diagonal(x = f$area * scale[root])
    d <- nmf_profile_change(amplitudes(a, old$scale), amplitudes(b, new$scale))
    change[at] <- pmax(change[at], d)
  }
  change
}

# Loss of deleting one occurrence with all other compounds in the offset.
# Read only occupied scans; observations, weights and offsets stay on their
# original scale. At each shared root we can remove individual occurrences
# or retire the entire normalized spectrum and save its m-1 degrees of freedom.
sinmf_delete_proposal <- function(state, price) {
  fits <- lapply(state$paths, readRDS)
  roots <- lapply(fits, function(f) state$map[f$gid])
  ids <- sort(unique(unlist(roots)))
  area <- numeric(nrow(state$P))
  for (s in seq_along(fits)) area[roots[[s]]] <- area[roots[[s]]] + fits[[s]]$area * state$scale[roots[[s]]]
  for (id in ids[order(area[ids])]) {
    k <- vapply(roots, function(g) match(id, g, nomatch = 0L), 0L)
    samples <- which(k > 0L)
    loss <- vapply(samples, function(s) {
      f <- fits[[s]]
      at <- which(f$H[, k[s]] > 0)
      if (!length(at)) return(0)
      y <- sinmf_sample_data(state$store, s, rows = at)
      P <- state$P[roots[[s]],,drop=FALSE] * (f$area * state$scale[roots[[s]]])
      exposure <- rep_len(f$exposure, nrow(f$H))[at]
      mu <- as.matrix(f$H[at,,drop=FALSE] %*% P) * exposure +
        (f$background[at,,drop=FALSE] %*% f$alpha) * rep_len(f$raw_exposure, nrow(f$H))[at] +
        sinmf_observation_block(f$offset, at)
      off <- mu - (as.numeric(f$H[at,k[s]]) * exposure) %o% P[k[s],]
      f$weights <- sinmf_observation_block(f$weights, at)
      2 * (sinmf_loglik(f, y, mu) - sinmf_loglik(f, y, pmax(off, 0)))
    }, 0)
    if (any(!is.finite(loss))) next
    remove <- if (sum(pmax(loss - price, 0)) < price * (ncol(state$P) - 1L))
      samples else samples[loss < price]
    for (s in remove) fits[[s]]$area[k[s]] <- 0
  }
  alive <- lapply(fits, function(f) which(f$area > 0))
  if (sum(lengths(alive)) == state$L) return(state)
  live <- sort(unique(unlist(Map(function(g, j) g[j], roots, alive))))
  map <- integer()
  for (s in seq_along(fits)) {
    j <- alive[[s]]
    f <- fits[[s]]
    f$H <- f$H[, j, drop = FALSE]
    f$area <- f$area[j]
    f$gid <- length(map) + seq_along(j)
    map <- c(map, match(roots[[s]][j], live))
    if (!is.null(f$profile_limits)) f$profile_limits <- f$profile_limits[,j,drop=FALSE]
    state <- sinmf_checkpoint(state, s, f, "support")
  }
  state$map <- map
  state$P <- state$P[live,,drop=FALSE]
  state$scale <- state$scale[live]
  state$L <- sum(lengths(alive))
  state
}

sinmf_end_refine <- function(state, control, q) {
  ids <- sort(unique(state$map))
  before <- sinmf_score(state)
  step <- function(f,j,i) sinmf_optimize(f, cycles = 1L, roots = j, q = q)
  run <- nmf_sweep_schedule(state, ids, control$refinement_sweeps, control$selective, step,
                             identity, function(a,b) sinmf_factor_changes(a,b,ids))
  candidate <- run$state
  best <- sinmf_score(candidate)
  if (is.finite(best) && best <= before) state <- candidate else best <- before
  retained <- if (control$refinement_sweeps && best < before) "ordinary sweeps" else "input"
  if (control$penalty_search && state$lambda > 0 && length(ids)) {
    candidate <- sinmf_delete_proposal(state, 4 * state$lambda)
    if (candidate$L < state$L) {
      candidate <- sinmf_optimize(candidate, cycles = 1L, q = q, record = FALSE)
      score <- sinmf_score(candidate)
      if (is.finite(score) && score < best) {
        candidate$history <- rbind(state$history, data.frame(stage = "support search", cycle = 1L,
          objective = score, delta = score - best, seconds = NA_real_))
        state <- candidate
        retained <- "fourfold deletion"
        best <- score
      }
    }
  }
  state$refinement <- list(before = before, after = best, retained = retained,
                            visited = run$visited, possible = run$possible)
  state
}
