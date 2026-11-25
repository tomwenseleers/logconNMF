# An alternate PSF supplies independent proposals from the same observations.
# A supplied design may come from calibration standards. TRUE estimates a
# design from the observed sample using the package's blind PSF estimator.
nmf_restart_basis <- function(x, y, time, family, ncores) {
  if (isTRUE(x))
    x <- estimatePSF(pmax(as.matrix(y),0), time=time, family=family,
                     nwindows=1L, precluster=1L, ncores=ncores)$X
  if (length(dim(x)) != 2L || nrow(x) != nrow(y) || ncol(x) < 1L)
    stop("alternative_X must have one row per original scan and at least one column.")
  x
}

# Pool integrated signal on common time bins, preserving each sample's total.
# The shared PSF is transferred to each native grid before fitting spectra.
# Only one raw sample and one vector of scan totals per sample are needed.
nmf_pooled_restart <- function(samples, family, dispersion, correlation, lambda, ncores) {
  counts <- times <- vector("list",length(samples))
  phi <- numeric(length(samples))
  for (i in seq_along(samples)) {
    s <- samples[[i]]
    y <- if (is.character(s$y)) readRDS(s$y) else s$y
    times[[i]] <- as.numeric(s$time)
    if (length(times[[i]]) != NROW(y) || NROW(y) < 9L ||
        any(!is.finite(times[[i]])) || any(diff(times[[i]]) <= 0)) stop("Invalid sample data/time.")
    counts[[i]] <- as.numeric(Matrix::rowSums(y))
    spec <- s$noise_spec
    if (!is.null(s$fit)) {
      f <- if (is.character(s$fit)) readRDS(s$fit) else s$fit
      if (!is.null(f$model)) f <- f$model
      spec <- f$spec
    }
    if (is.null(spec)) {
      d <- if (is.null(s$dispersion)) dispersion else s$dispersion
      spec <- estimateSignalNoise(y,family,dispersion=d,correlation=correlation,ncores=ncores)$spec
      samples[[i]]$noise_spec <- spec
    }
    phi[i] <- spec$dispersion
  }
  n <- as.integer(round(stats::median(lengths(times))))
  edges <- lapply(times,sinmf_edges)
  bins <- seq(min(vapply(edges,min,0)),max(vapply(edges,max,0)),length.out=n+1L)
  t <- utils::head(bins,-1L)+diff(bins)/2
  total <- Reduce(`+`,Map(function(e,v) diff(stats::approx(e,c(0,cumsum(v)),bins,rule=2)$y), edges,counts))
  weight <- vapply(counts,sum,0)
  if (!any(weight > 0)) return(lapply(samples,function(s) {s$alternate_design <- NULL; s}))
  common <- estimatePSF(matrix(total,n,1L),time=t,family=stats::quasipoisson("identity"),method="blind",
    dispersion=stats::weighted.mean(phi,weight),
    precluster=FALSE,nwindows=2L,ncores=ncores,
    control=list(maxiter=3L,fit=list(IC=if (is.character(lambda)) lambda else "gic")))
  for (i in seq_along(samples))
    samples[[i]]$alternate_design <- estimatePSF(matrix(counts[[i]],ncol=1L),time=times[[i]],
      method="transfer",standard=common,ncores=ncores)$X
  samples
}

# Preserve the model, noise parameters and penalty of the primary fit when
# comparing a general-family alternate. Only the PSF initialization changes.
nmf_alternate_single <- function(X, y, fit, family, lambda, precluster, weights, offset, control, ncores) {
  candidate <- splogconNMF(X, y, family=family, lambda=lambda, background=fit$background,
    dispersion=fit$spec$dispersion, theta=if (fit$spec$kind == 3L) fit$spec$theta,
    precluster=precluster, weights=weights, offset=offset, correlation="none",
    control=control, ncores=ncores)
  candidate$spec <- fit$spec
  candidate
}

slcn_restart <- function(fit, challenger, y, refinement) {
  start <- proc.time()[["elapsed"]]
  alternate <- fit
  alternate$H <- challenger$H
  alternate$S <- challenger$S
  alternate$alpha <- challenger$alpha
  alternate$state_ids <- alternate$source_ids <- NULL
  alternate <- slcn_sync(alternate,y)
  if (refinement$refinement_sweeps || refinement$penalty_search)
    alternate <- slcn_end_refine(alternate,y,refinement)
  score <- alternate$objective
  before <- fit$objective
  if (score < before) {
    alternate$history <- rbind(fit$history, data.frame(stage="alternative PSF", seconds=proc.time()[["elapsed"]]-start,
      delta=score-before, accepted=TRUE, compounds=alternate$rank))
    alternate$complete_delta <- score-fit$initial_objective
    alternate$ledger_error <- alternate$complete_delta-sum(alternate$history$delta[alternate$history$accepted])
    fit <- slcn_temporal_order(alternate)
  }
  fit$restart <- list(before=before, alternative=score, after=fit$objective)
  fit
}
