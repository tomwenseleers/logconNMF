# Both blind fits have their own registration. Heights live on the same
# original scans. Only the search bands need converting to recipient coordinates.
cc_restart_bands <- function(bands, old_warps, data) {
  if (!nrow(bands)) return(bands)
  mapped <- Map(function(w, x) matrix(stats::approx(w, x$warp, xout = as.numeric(bands),
    rule = 2)$y, ncol = 2L), old_warps, data)
  cbind(Reduce(pmin, lapply(mapped, function(b) b[,1L])),
        Reduce(pmax, lapply(mapped, function(b) b[,2L])))
}

cc_restart_search <- function(z, donor) {
  score <- function() cc_dimension(z$data, z$p, cc_D(z$data), z$price)$objective
  before <- score()
  incumbent <- cc_checkpoint(z)
  cc_restore(z, donor)
  alternate <- score()
  if (alternate < before) {
    donor <- incumbent
  } else cc_restore(z, incumbent)
  initial <- min(before, alternate)
  primary <- nrow(z$p)
  p <- rbind(z$p, donor$p)
  bands <- rbind(z$bands, donor$bands)
  h <- Map(function(x, other) cbind(as(x$h, "dgCMatrix"), other), z$data, donor$h)
  active <- seq_len(nrow(p)) <= primary
  group <- seq_along(active)
  df <- rowSums(p > 0) - 1L
  for (s in seq_along(h)) {
    hs <- as.matrix(h[[s]])
    df <- df + cc_faces(z$data[[s]]$t, hs, TRUE)
    for (ix in cc_overlap_groups(hs)) {
      join <- unique(group[ix])
      if (length(join) > 1L) group[group %in% join] <- min(join)
    }
  }
  groups <- split(seq_along(group), group)
  swaps <- 0L
  repeat {
    changed <- FALSE
    for (ix in groups) {
      from <- ix[active[ix]]
      to <- ix[!active[ix]]
      ions <- which(colSums(p[ix,,drop=FALSE]) > 0)
      delta <- z$price * (sum(df[to]) - sum(df[from]))
      updates <- vector("list", length(h))
      for (s in seq_along(h)) {
        x <- z$data[[s]]
        rows <- which(Matrix::rowSums(h[[s]][,ix,drop=FALSE]) > 0)
        if (!length(rows)) next
        component <- function(j) as.matrix(h[[s]][rows,j,drop=FALSE] %*% p[j,ions,drop=FALSE])
        mu <- x$mu[rows,ions,drop=FALSE]
        next_mu <- pmax(mu - component(from) + component(to), 1e-8)
        delta <- delta + cc_diff(x$y[rows,ions,drop=FALSE], next_mu, mu, x$phi)
        updates[[s]] <- list(rows = rows, mu = next_mu)
      }
      if (is.finite(delta) && delta < -1e-5) {
        for (s in which(lengths(updates) > 0L))
          cc_set_mu(z$data[[s]], updates[[s]]$rows, ions, updates[[s]]$mu)
        active[ix] <- !active[ix]
        swaps <- swaps + 1L
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  z$p <- p[active,,drop=FALSE]
  z$bands <- bands[active,,drop=FALSE]
  for (s in seq_along(h)) {
    z$data[[s]]$h <- as.matrix(h[[s]][,active,drop=FALSE])
    cc_means(z$data[[s]], z$p)
  }
  after <- score()
  if (after > initial + 1e-5) stop("Local PSF exchange increased the original-channel objective.")
  if (swaps) {
    saved <- cc_checkpoint(z)
    for (i in seq_len(30L)) cc_end_step(z, seq_len(nrow(z$p)), i)
    cc_support_sweep(z$data,z$p,order(colSums(cc_mass(z$data))),z$price)
    cc_trim(z)
    if (score() > after) cc_restore(z, saved)
  }
  if (z$trace) cc_history(z, "alternative PSF", 1L, "local support search", cc_D(z$data))
  list(before = before, alternative = alternate, swaps = swaps, after = score())
}
