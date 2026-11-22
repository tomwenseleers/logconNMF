#' Fit a log-concave factorization with automatic rank selection
#'
#' Estimate ordered, nonnegative log-concave profiles and channel loadings from
#' one signal or several related samples. The information criterion selects the
#' number of components. Simultaneous fits share spectra while retaining the
#' original scan grids and estimating a separate profile and area for each
#' occurrence.
#'
#' @name logconNMF
#' @order 1
#' @param y A numeric scans-by-channels matrix or vector for one sample, or a
#'   list of at least two sample descriptions for a simultaneous fit. Each
#'   description contains `y` and increasing `time`, with optional `X`,
#'   `background`, `fit`, `weights`, `offset`, `exposure`, `dispersion`,
#'   `theta`, `sample_id`, `group`, `small_width` and `edges`. Supply
#'   observations, weights, offsets or initial fits as RDS paths to
#'   read them from disk. Samples share channel names and order. List names
#'   supply sample IDs when the descriptions omit `sample_id`.
#' @param X Candidate profile matrix for a single sample, such as the sparse
#'   design from [estimatePSF()]. In a simultaneous fit, supply it in each sample
#'   description or let the simultaneous initializer estimate it.
#' @param lambda An information criterion (`"gic"`, `"aic"`, `"bic"`,
#'   `"ebic"` or `"mbic"`) or a nonnegative numeric penalty per selection
#'   degree of freedom, before any temporal information adjustment. With the general solver and
#'   \eqn{m} channels, a single component incurs \eqn{m\lambda}; in a
#'   simultaneous fit a shared spectrum incurs \eqn{(m-1)\lambda} and an
#'   occurrence incurs \eqn{\lambda}.
#' @param precluster Maximum channel groups or a saved [channelGroups()] map.
#'   `"auto"` uses the original channels for a single sample and up to 32
#'   common groups for a simultaneous fit. `FALSE` or `NULL` disables aggregation.
#'   Final fits use every original channel.
#' @param correlation Lag-one working information adjustment. `NULL` uses the
#'   fitter's default (`"auto"` for one sample and `"none"` for simultaneous fits),
#'   or `control$correlation` if supplied for a simultaneous fit.
#' @param control Named controls. See [splogconNMF()] for one matrix. The conditional count solver uses up to
#'   12 initial sweeps, stopping at a relative improvement of 0.001. For a
#'   general simultaneous fit, `cycles=10L` and `initial_cycles=5L` set refinement limits,
#'   `rounds=3L` sets the number of grouping rounds, and `q=64L` sets the
#'   initial profile update resolution. After a cross-sample match, each
#'   local profile is refined on all original channels using the hard shape
#'   constraint of [fitLogConcaveProfile()]. Its window extends two local
#'   PSF full widths at half maximum on either side of the initial support.
#'   The default `rt_tolerance=NULL` and
#'   `alignment_window=NULL` derive matching windows from the narrow-peak
#'   widths. `ridge=TRUE` checks the profile ridge across samples. The
#'   remaining defaults are `censor=0.5`, `correlation="none"`,
#'   `workers=NULL` and `cache=TRUE`. `initial` passes controls to
#'   the independent initial fits and defaults to
#'   `list(initial_relaxed=FALSE)`; single-sample fits default to a relaxed
#'   grouped-support initialization. `workers` limits concurrent samples
#'   within `ncores`. Conditional grouped regressions use a memory-aware worker cap.
#'   For the general solver the cache retains aggregated observations and at most
#'   one original count matrix up to 64 MiB plus its likelihood normalizers per
#'   worker. `cache=FALSE` disables persistent input caches and temporary
#'   whole-sample profile caches, retaining local-window disk reads.
#'   For memory-limited series, `workers=1L, cache=FALSE` avoids concurrent
#'   sample workspaces and retained observation caches. Set
#'   `initial=list(initial_solver=list(solver="eigen", method=3L, precond=1L))`
#'   to request iterative component selection for every sample. The solver's
#'   `max.memory` budget covers Gram/factor workspace, not total process RAM.
#' @param ncores Total core budget. `NULL` uses physical cores for one sample
#'   and logical cores for simultaneous fits. Grouped PSF regressions divide this
#'   budget between samples and native threads. The conditional count solver then
#'   reuses it for sample-specific profile and background updates and disjoint
#'   ion columns of the shared spectra, with one BLAS thread per worker. Set
#'   `ncores=1L` for a serial fit.
#' @param path New directory for the sample store and fitted factors in a
#'   simultaneous fit. `NULL` uses a temporary directory.
#' @param align Estimate monotone time warps for a simultaneous fit.
#' @param method Spectral distance for matching simultaneous samples, either
#'   `"cosine"` or `"hellinger"`.
#' @param object A fitted model returned by `logconNMF()`.
#'   `sampleNMF()` takes a simultaneous model.
#' @param x A fitted simultaneous model for plotting.
#' @param sample Sample index or ID. For a single matrix use `1L` or
#'   `"sample1"`. In `coef()`, `NULL` returns the common spectra for a
#'   simultaneous fit. In `plot()`, `NULL` displays all samples.
#' @param normalize `"none"` returns amplitude-scaled loadings. `"l1"`
#'   returns unit-sum rows.
#' @param type Plot sample areas, common spectra or aligned compound profiles.
#' @param compounds Compound indices to display in a simultaneous plot.
#' @param ... Further arguments passed to the underlying extractor.
#' @inheritParams splogconNMF
#' @aliases simultaneousNMF
#' @details
#' **Rank-selection objective.** The fit approximately minimizes
#' \deqn{-2\sum_s \ell_{\nu_s}(M_s;Y_s) + \kappa\lambda d.}
#' The weighted log-likelihood \eqn{\ell_{\nu_s}} uses the fitted mean
#' \eqn{M_s} on the original channels of sample \eqn{s}. For quasipoisson
#' it is the Poisson log score divided by the sample's estimated or supplied
#' dispersion. The optional working temporal adjustment is \eqn{\kappa},
#' equal to one under independence. A numeric `lambda` sets the penalty per
#' selection degree of freedom. The GIC preset uses \eqn{\log(p)\log\log(n)},
#' with the PSF dictionary width and number of scans for a single matrix.
#' For a sample list the preset uses \eqn{p=n=\sum_s n_s}.
#'
#' For the conditional count solver the complexity is the sum of the number of
#' positive ions minus one per normalized spectrum and the dimension of every
#' positive log-profile face (two endpoints plus interior slope changes).
#' \code{fit$conditional_objective} reports quasi-deviance and both dimensions.
#' An absent occurrence has zero cost. Finite support proposals and alternating
#' conditional updates give a local optimum candidate. For the general solver,
#' for one matrix with \eqn{m} channels and \eqn{K} components,
#' \eqn{d=mK}. For a simultaneous model, \eqn{d=(m-1)G+L}, where \eqn{G} is
#' the number of shared unit-sum spectra and \eqn{L} counts the detected
#' sample–component occurrences. Each spectrum contributes \eqn{m-1}
#' coordinates and each occurrence contributes one amplitude. Setting
#' \eqn{G=L=K} gives the single-sample count. Profiles have contiguous
#' support and a concave logarithm on the positive interval. Profile shapes
#' and smooth backgrounds carry no additional rank charge.
#'
#' The simultaneous `objective` and reported information criteria such as [GIC()]
#' also count estimated noise parameters. The single-sample `objective`
#' drops that constant charge for a fixed noise specification. Local
#' refinement and fusion, deletion or addition of candidate components
#' provide the approximate rank search.
#' Final single-sample refinement also tries broader local starts for nearly
#' coincident profiles. Compatible Poisson/quasipoisson fits use their existing
#' channel groups for preliminary screening; completed candidates and controls
#' are always refined and checked on every original channel. Other cases retain
#' full-channel screening. This is enabled
#' by default. For one matrix, `control=list(profile_restarts=FALSE)` disables
#' the pass; simultaneous fits pass this option through `control$initial`.
#'
#' For a simultaneous fit, put `background`, `weights`, `offset` and any initial
#' `fit` in the sample descriptions. `path = NULL` uses a new temporary
#' sample store. Supply a persistent path when saving a simultaneous fit for later
#' use, and retain the store with the fitted object.
#'
#' Numeric penalties in `logconNMF()` have the same scale in both modes.
#' To reproduce a single-sample [splogconNMF()] call using a numeric compound
#' penalty, divide that penalty by the number of channels.
#'
#' @section Simultaneous models:
#' Ignoring backgrounds, exposures and offsets, the sample means have the
#' shared-factor form \eqn{M_s = Q_s P}, with one common spectral factor
#' \eqn{P} and sample-specific profile/area factors \eqn{Q_s}. This follows the
#' simultaneous-NMF terminology of
#' [Cichocki et al. (2009)](https://doi.org/10.1002/9780470747278) and
#' [Lee et al. (2012)](https://doi.org/10.1371/journal.pone.0048238).
#'
#' Samples can supply separate background bases, known exposures and noise
#' parameters. With a time warp \eqn{u=w_s(t)}, binwise Jacobian \eqn{J},
#' and exposure \eqn{e}, the fitted mean on the original grid is
#' \deqn{\mu_{sij}=e_{si}\{J_{si}\sum_g a_{sg}h_{sg}(u_{si})p_{gj}
#' +(B_sA_s)_{ij}\}+o_{sij}.}
#' Spectra sum to one and the profiles integrate to one on the aligned grid.
#' Unit exposure models the supplied intensities. Use known integration
#' durations as exposures for integrated counts. Scan spacing and
#' `small_width` use the same time units.
#'
#' Reciprocal spectral matches provide leading-edge alignment anchors.
#' Compatible occurrences share spectra after checking overlap and profile
#' shape in area order. At most one occurrence of a compound enters each
#' sample. Local profile, spectral and background updates compare the
#' full original-channel likelihood. Undetected occurrences have `NA`
#' areas. Noise scales are held fixed when comparing groupings. A negative
#' NB2 `theta` allows bounded, likelihood-improving shape updates.
#'
#' The common channel map is learned from samples aggregated on common time
#' bins. The conditional solver holds original counts and means in memory for
#' joint updates and releases the initial PSF workers before this stage.
#' the general solver retains only the aggregated observations in its shared
#' cache. Spectral, background and score passes stream original-channel blocks.
#' Poisson/quasipoisson workers may additionally retain one original sample
#' up to 64 MiB and an equally sized log-factorial cache. Larger samples use
#' local reads. These bounded caches are released with the runtime, including
#' error exits; they do not grow with the number of samples. Absolute likelihoods
#' retain their normalizers and local acceptance checks retain all original ions.
#' Other profile passes load one assigned sample per active worker.
#' Keep the store with the saved fitted object for reconstruction and plots.
#'
#' @section Extract one sample:
#' `sampleNMF(object, sample)` returns a single-sample factorization view.
#' The view contains unit-sum original-grid profiles in `H`, matching
#' amplitude-scaled spectra in `S`, its fitted background, integrated areas
#' in `integrated_area`, and consensus IDs in `compound_id`. Alignment and
#' exposure are included in the profiles and spectra. Use `coef()`,
#' `fitted()` and `GIC()` on the returned view. The view carries no observed
#' matrix.
#'
#' @return A fitted object with estimated `rank`, the number of single-sample
#'   components or shared components. Use [coef()], [fitted()] and [GIC()] in
#'   either mode. `fitted()` returns one sample on its original grid. A
#'   simultaneous fit includes consensus spectra `P`, sample-by-compound `areas`,
#'   `compounds`, `alignment`, `paths`, `store`, `channel_groups` and
#'   `parallel` thread allocation.
#' @examples
#' X <- convMatrix(convKernel("gaussian", 9L), 50L)
#' signal <- cbind(700 * X[, 16L] + 350 * X[, 35L],
#'                 150 * X[, 16L] + 500 * X[, 35L])
#' set.seed(12)
#' Y <- matrix(rpois(length(signal), signal), nrow(signal))
#' fit <- logconNMF(Y, X = X, ncores = 1L)
#' fit$rank
#' coef(fit, normalize = "l1")
#'
#' second <- matrix(rpois(length(signal), 1.4 * signal), nrow(signal))
#' samples <- list(first = list(y = Y, time = seq_len(nrow(Y)), X = X),
#'   second = list(y = second, time = seq_len(nrow(Y)), X = X))
#' joint <- logconNMF(samples, ncores = 1L)
#' joint$rank
#' joint$areas
#' fitted(joint, sample = "second")[1:3, ]
#' second_fit <- sampleNMF(joint, sample = "second")
#' second_fit$compound_id
#' coef(second_fit)
#' @md
#' @param algorithm Solver for count data. `"conditional"` alternates hard
#'   log-concave profiles and sparse spectra on the original ion counts.
#'   Multisample fits can reassign mixed spectra to recurring components. A local
#'   refit must improve the same original-channel information criterion.
#'   `"classic"` uses the general identity-link solver. `"auto"` uses the
#'   built-in conditional solver for supported inputs.
#' @param time Increasing scan coordinates for a single sample. Defaults to scan
#'   index. Supply native time in each descriptor for a simultaneous fit.
#' @param objective.trace Record the penalized quasi-deviance after each
#'   conditional block. The default avoids the cost of those measurements.
#' @export
logconNMF <- function(y, X = NULL, family = poisson("identity"), lambda = "gic",
                      background = NULL, dispersion = NULL, theta = NULL,
                      precluster = "auto", weights = NULL, offset = NULL,
                      correlation = NULL, ncores = NULL, control = list(),
                      start = NULL, path = NULL, align = TRUE,
                      method = c("cosine", "hellinger"), verbose = FALSE,
                      algorithm = c("auto", "conditional", "classic"), time = NULL,
                      objective.trace = FALSE) {
  algorithm <- match.arg(algorithm)
  if (!is.logical(objective.trace) || length(objective.trace) != 1L || is.na(objective.trace))
    stop("objective.trace must be TRUE or FALSE.")
  begun <- proc.time()[["elapsed"]]
  if (is.numeric(lambda) && (length(lambda) != 1L || is.complex(lambda) ||
      !is.finite(lambda) || lambda < 0)) {
    stop("Use one finite, nonnegative penalty or an information criterion.", call. = FALSE)
  }
  if (!is.list(control)) stop("Use a named control list.", call. = FALSE)
  simultaneous <- is.list(y) && !is.data.frame(y)
  if (identical(precluster, "auto")) precluster <- if (simultaneous) 32L else NULL
  if (simultaneous) {
    if (!is.null(time)) stop("Supply time in each sample description.")
    if (length(y) < 2L) stop("Supply a matrix for one sample or at least two sample descriptions.", call. = FALSE)
    if (!is.null(X) || !is.null(background) || !is.null(weights) ||
        !is.null(offset) || !is.null(start)) {
      stop("Supply X, background, weights, offset and initial fits in the sample descriptions.",
        call. = FALSE)
    }
    samples <- lapply(seq_along(y), function(i) {
      sample <- y[[i]]
      if (!is.list(sample) || is.null(sample$y) || is.null(sample$time)) {
        stop("Each sample description needs y and time.", call. = FALSE)
      }
      name <- names(y)[i]
      if (is.null(sample$sample_id) && length(name) && !is.na(name) && nzchar(name)) {
        sample$sample_id <- name
      }
      sample
    })
    if (!is.null(correlation)) control$correlation <- correlation
    if (is.null(path)) path <- tempfile("logconNMF_")
    classic_controls <- intersect(names(control), c("rounds", "initial_cycles", "ridge", "q"))
    if (length(classic_controls) && algorithm == "conditional")
      stop("These controls are specific to the general solver: ", paste(classic_controls, collapse = ", "))
    conditional <- !length(classic_controls) && conditional_route(algorithm, family, dispersion, theta, samples)
    fit <- simultaneousNMF(samples, family = family, lambda = lambda, path = path,
      precluster = precluster, dispersion = dispersion, theta = theta,
      method = match.arg(method), align = align, ncores = ncores, control = control,
      verbose = verbose, conditional = conditional, trace = objective.trace)
  } else {
    if (is.null(X)) stop("Supply X for a single-sample fit.", call. = FALSE)
    if (!is.null(path) || !missing(align) || !missing(method)) {
      stop("path, align and method apply to simultaneous fits.", call. = FALSE)
    }
    scalar <- lambda
    conditional <- conditional_route(algorithm, family, dispersion, theta, list(list(weights=weights, offset=offset)))
    if (is.numeric(lambda)) lambda <- lambda * NCOL(y)
    if (is.null(correlation)) correlation <- "auto"
    fit <- splogconNMF(X, y, family = family, background = background, lambda = lambda,
      dispersion = dispersion, theta = theta, precluster = precluster, weights = weights,
      offset = if (conditional && is.null(background)) 1e-8 else offset,
      correlation = correlation, ncores = ncores, control = control,
      start = start, verbose = verbose)
    if (conditional) {
      price <- fit$spec$information_factor * if (is.numeric(scalar)) scalar else
        presetLambda(scalar, p = fit$p, n = fit$n)
      fit <- conditional_single(fit, as.matrix(y), time, price, objective.trace)
      fit$seconds["total"] <- proc.time()[["elapsed"]] - begun
    }
  }
  if (is.null(fit$algorithm)) fit$algorithm <- "classic"
  class(fit) <- c("logconNMF_fit", class(fit))
  fit
}

logconNMF_check_sample <- function(sample) {
  if (!is.null(sample) && !isTRUE(sample %in% c("1", "sample1"))) {
    stop("Select sample 1 or 'sample1' for a single-sample fit.", call. = FALSE)
  }
}

#' @rdname logconNMF
#' @export
coef.logconNMF_fit <- function(object, sample = NULL, normalize = c("none", "l1"), ...) {
  normalize <- match.arg(normalize)
  if (inherits(object, "simultaneousNMF")) {
    return(coef.simultaneousNMF(object, sample = sample, normalize = normalize, ...))
  }
  logconNMF_check_sample(sample)
  loadings <- coef.splogconNMF(object, ...)
  if (normalize == "l1" && nrow(loadings)) {
    totals <- rowSums(loadings)
    totals[totals == 0] <- 1
    loadings <- loadings / totals
  }
  loadings
}

#' @rdname logconNMF
#' @export
fitted.logconNMF_fit <- function(object, sample = 1L, ...) {
  if (inherits(object, "simultaneousNMF")) return(fitted.simultaneousNMF(object, sample = sample, ...))
  logconNMF_check_sample(sample)
  fitted.splogconNMF(object, ...)
}
