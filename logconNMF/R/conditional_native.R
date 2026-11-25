conditional_native <- function(parts, group, price, trace, cycles, exchanges, ncores, reference=NULL,
                               refinement=nmf_refinement_control(list()), alternate=NULL, checkpoint=FALSE) {
  if (!length(parts) || !length(group) || anyNA(group) || any(group < 1) ||
      any(group != as.integer(group)) || !is.finite(price) || price < 0)
    stop("Invalid channel groups, penalty or empty count-factorization input.")
  z <- new.env(parent = emptyenv())
  z$data <- if (is.null(reference)) lapply(parts, cc_sample, group = group, price = price, threads = ncores) else
    Map(function(part, ref) cc_reference_sample(part,ref,group,price,ncores),parts,reference)
  z$price <- price
  z$trace <- isTRUE(trace)
  z$threads <- as.integer(ncores)
  z$history <- list()
  dictionary <- if (is.null(reference)) cc_dictionary(z$data, price, ncores) else cc_reference_dictionary(z$data)
  z$p <- dictionary$p
  z$bands <- dictionary$bands
  for (x in z$data) {
    rm(list = c("init", "initial_p", "norms"), envir = x)
    cc_means(x, z$p)
  }
  selective <- refinement$selective_initial
  initial_work <- cc_optimize(z, cycles, "factor updates", TRUE, selective)
  cc_trim(z)
  if (nrow(z$p)) {
    for (i in seq_len(exchanges)) {
      before <- cc_dimension(z$data, z$p, cc_D(z$data), price)$objective
      moves <- cc_exchange(z)
      D <- cc_D(z$data)
      after <- cc_dimension(z$data, z$p, D, price)$objective
      if (after > before + 1e-04)
        stop("Support exchange increased the joint original-ion criterion.")
      cc_history(z, "support exchange", i, "exchange", D)
      if (!moves)
        break
      cc_optimize(z, 10L, paste("conditional refit", i), selective = selective)
    }
    previous <- cc_dimension(z$data, z$p, cc_D(z$data), price)$objective
    change <- cc_shared_exchange(z)
    if (length(z$data) > 1L) {
      change <- change + cc_split_spectra(z)
      change <- change + cc_shared_exchange(z, TRUE)
      change <- change + cc_refit_overlaps(z)
    }
  } else {
    previous <- cc_dimension(z$data, z$p, cc_D(z$data), price)$objective
    change <- 0
  }
  cc_trim(z)
  for (x in z$data) cc_means(x, z$p)
  D <- cc_D(z$data)
  objective <- cc_dimension(z$data, z$p, D, price)
  if (abs(objective$objective - (previous + change)) > 1e-04 + 1e-08 * abs(objective$objective))
    stop("Local updates disagree with the complete original-ion objective.")
  cc_history(z, "local joint refit", 1L, "overlap", D)
  report <- NULL
  if (refinement$refinement_sweeps || refinement$penalty_search) {
    report <- if (refinement$penalty_search && refinement$refinement_sweeps == 0L &&
                  is.null(alternate) && !checkpoint) {
      cc_local_search(z, objective$objective, refinement$selective)
    } else cc_end_refine(z, refinement)
    report$initial <- initial_work
    cc_trim(z)
    objective <- cc_dimension(z$data, z$p, cc_D(z$data), price)
  }
  if (checkpoint) {
    snap <- cc_checkpoint(z)
    snap$warps <- lapply(z$data, function(x) x$warp)
    return(snap)
  }
  restart <- NULL
  if (!is.null(alternate)) {
    alternate$bands <- cc_restart_bands(alternate$bands, alternate$warps, z$data)
    restart <- cc_restart_search(z, alternate, selective)
    objective <- cc_dimension(z$data, z$p, cc_D(z$data), price)
  }
  idx <- order(rowMeans(z$bands))
  factors <- lapply(z$data, function(x) {
    h <- x$h[, idx, drop = FALSE]
    list(H = as(h, "dgCMatrix"), alpha = x$alpha, shape_df = cc_faces(x$t, h, TRUE))
  })
  list(P = z$p[idx, , drop = FALSE], factors = factors, objective = objective, refinement = report, restart = restart, trace = if (z$trace) do.call(rbind,
    z$history) else NULL)
}
