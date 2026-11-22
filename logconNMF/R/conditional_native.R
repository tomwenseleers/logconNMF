conditional_native <- function(parts, group, price, trace, cycles, exchanges, ncores) {
  if (!length(parts) || !length(group) || anyNA(group) || any(group < 1) ||
      any(group != as.integer(group)) || !is.finite(price) || price < 0)
    stop("Invalid channel groups, penalty or empty count-factorization input.")
  z <- new.env(parent = emptyenv())
  z$data <- lapply(parts, cc_sample, group = group, price = price, threads = ncores)
  z$price <- price
  z$trace <- isTRUE(trace)
  z$threads <- as.integer(ncores)
  z$history <- list()
  dictionary <- cc_dictionary(z$data, price, ncores)
  z$p <- dictionary$p
  z$bands <- dictionary$bands
  for (x in z$data) {
    rm(list = c("init", "initial_p", "norms"), envir = x)
    cc_means(x, z$p)
  }
  cc_optimize(z, cycles, "factor updates", TRUE)
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
      cc_optimize(z, 10L, paste("conditional refit", i))
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
  idx <- order(rowMeans(z$bands))
  factors <- lapply(z$data, function(x) {
    h <- x$h[, idx, drop = FALSE]
    list(H = as(h, "dgCMatrix"), alpha = x$alpha, shape_df = cc_faces(x$t, h, TRUE))
  })
  list(P = z$p[idx, , drop = FALSE], factors = factors, objective = objective, trace = if (z$trace) do.call(rbind,
    z$history) else NULL)
}
