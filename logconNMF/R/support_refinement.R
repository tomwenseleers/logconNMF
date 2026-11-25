# The temporary support price proposes a model. Only the original, ungrouped
# observations and the original price decide whether to retain it.
nmf_refinement_control <- function(control) {
  defaults <- list(refinement_sweeps = 0L, penalty_search = FALSE, selective = TRUE)
  out <- utils::modifyList(defaults, control[intersect(names(defaults), names(control))], keep.null = TRUE)
  n <- out$refinement_sweeps
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 0 ||
      n > .Machine$integer.max || n != floor(n)) stop("refinement_sweeps must be a nonnegative integer.")
  for (key in c("penalty_search", "selective"))
    if (!is.logical(out[[key]]) || length(out[[key]]) != 1L || is.na(out[[key]]))
      stop(key, " must be TRUE or FALSE.")
  out$refinement_sweeps <- as.integer(n)
  out
}

# At regular full sweeps check both the unit spectrum and each sample's
# amplitude-scaled profile. An intervening change in a neighbour or the
# background can therefore bring a previously stationary factor back.
nmf_sweep_schedule <- function(state, ids, n, selective, step, snapshot, change) {
  use <- ids
  visited <- 0
  for (i in seq_len(n)) {
    audit <- selective && n >= 10L && (i == 1L || i %% 10L == 0L || i == n)
    if (audit) use <- ids
    if (audit && i < n) before <- snapshot(state)
    state <- step(state, use, i)
    visited <- visited + length(use)
    if (audit && i < n) use <- ids[change(before, state) >= 1e-4]
  }
  list(state = state, visited = visited, possible = as.double(n) * length(ids))
}

nmf_profile_change <- function(old, new) {
  as.numeric(Matrix::colSums(abs(old - new)) /
    pmax(Matrix::colSums(old), Matrix::colSums(new), 1e-12))
}

cc_refinement_snapshot <- function(z) {
  list(p = z$p + 0, h = lapply(z$data, function(x) as(x$h, "dgCMatrix")))
}

cc_refinement_change <- function(before, z) {
  change <- rowSums(abs(before$p - z$p))
  for (s in seq_along(z$data))
    change <- pmax(change, nmf_profile_change(before$h[[s]], as(z$data[[s]]$h, "dgCMatrix")))
  change
}

cc_checkpoint <- function(z) {
  snap <- cc_refinement_snapshot(z)
  snap$bands <- z$bands + 0
  snap$alpha <- lapply(z$data, function(x) x$alpha + 0)
  snap
}

cc_restore <- function(z, snap) {
  z$p <- snap$p + 0
  z$bands <- snap$bands + 0
  for (s in seq_along(z$data)) {
    x <- z$data[[s]]
    x$h <- as.matrix(snap$h[[s]])
    x$alpha <- snap$alpha[[s]] + 0
    x$bg <- pmax(x$z %*% x$alpha + 1e-8, 1e-8)
    cc_means(x, z$p)
  }
}

cc_end_step <- function(z, ids, i, birth = FALSE, record = FALSE) {
  ordered <- function() {
    if (length(ids) < 2L) return(ids)
    mass <- if (length(ids) == nrow(z$p)) colSums(cc_mass(z$data))[ids] else
      Reduce(`+`, lapply(z$data, function(x) colSums(x$h[,ids,drop=FALSE])))
    ids[order(-mass)]
  }
  cc_profile_sweep(z$data, z$p, ordered(), z$bands, z$price, birth, z$threads)
  cc_spectrum_sweep(z$data, z$p, ordered(), z$price, TRUE, z$threads)
  if (record && z$trace) cc_history(z, "terminal refinement", i, "spectra and background", cc_D(z$data))
  z
}

cc_end_refine <- function(z, control) {
  price <- z$price
  objective <- function() cc_dimension(z$data, z$p, cc_D(z$data), price)$objective
  drop <- function() cc_support_sweep(z$data, z$p, order(colSums(cc_mass(z$data))), z$price)
  ids <- function() seq_len(nrow(z$p))
  before <- best <- objective()
  retained <- "input"
  if (control$refinement_sweeps) saved <- cc_checkpoint(z)
  nh <- length(z$history)
  run <- nmf_sweep_schedule(z, ids(), control$refinement_sweeps, control$selective,
    function(z, j, i) cc_end_step(z, j, i, record = TRUE), cc_refinement_snapshot, cc_refinement_change)
  now <- if (control$refinement_sweeps) objective() else before
  if (is.finite(now) && now <= best) {
    best <- now
    if (control$refinement_sweeps) retained <- "ordinary sweeps"
  } else {
    cc_restore(z, saved)
    z$history <- utils::head(z$history, nh)
  }
  if (control$penalty_search && price > 0 && nrow(z$p)) {
    saved <- cc_checkpoint(z)
    z$price <- 4 * price
    drop()
    cc_contract(z)
    for (i in seq_len(6L)) {
      cc_end_step(z, ids(), i)
      drop()
    }
    cc_trim(z)
    z$price <- price
    for (i in seq_len(20L)) cc_end_step(z, ids(), i)
    drop()
    now <- objective()
    if (is.finite(now) && now < best) {
      best <- now
      saved <- cc_checkpoint(z)
      retained <- "fourfold deletion"
    }
    cc_end_step(z, ids(), 1L, birth = TRUE)
    for (i in seq_len(10L)) cc_end_step(z, ids(), i)
    drop()
    now <- objective()
    if (is.finite(now) && now < best) {
      best <- now
      retained <- "fourfold deletion and re-entry"
    } else cc_restore(z, saved)
    cc_trim(z)
    if (z$trace) cc_history(z, "support search", 1L, "original price", cc_D(z$data))
  }
  list(before = before, after = best, retained = retained, visited = run$visited, possible = run$possible)
}
