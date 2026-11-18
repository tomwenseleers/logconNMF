#' Fit a log-concave rate profile to an ordered multichannel signal
#'
#' Estimates a single nonnegative profile given a fixed spectrum and the
#' contribution of the background and other components. With scan exposure
#' \eqn{e_t}, profile \eqn{h_t} and spectrum \eqn{s_j}, the fitted mean is
#' \eqn{O_{tj} + e_t h_t s_j}. Pass a local window to restrict its support.
#'
#' @param y Numeric matrix of observations, with ordered scans in rows and
#'   channels in columns, or a vector for one channel. Supply integer counts
#'   for Poisson and negative binomial, and encode censored Gamma values as zero.
#' @param spectrum Known nonnegative channel loadings in column order. A
#'   spectrum summing to one gives a profile in total-ion rate units.
#' @param offset Known nonnegative contribution from other components and
#'   background, in observation units. A scalar, scan vector or matrix like y.
#' @param family Identity-link noise model accepted by [logconNMF()].
#' @param time Strictly increasing scan coordinates, optionally irregular.
#'   Defaults to the row numbers.
#' @param exposure Positive scan exposures, scalar or one per row.
#' @param weights Nonnegative observation weights, scalar, scan vector or
#'   matrix like y. Zero weights omit observations from the profile update.
#' @param start Optional nonnegative initial rate profile. A second start is
#'   obtained from the data conditional on the supplied offset.
#' @param dispersion Fixed dispersion, or the scale for fixed-scale Gamma.
#'   Defaults to the parameter in the family object, otherwise one. For NB1,
#'   one gives the Poisson limit. Estimate noise from the larger sample
#'   with [estimateSignalNoise()] when appropriate.
#' @param theta Fixed positive size for NB2, supplied here or in the family
#'   object. The profile fit conditions on the noise parameters.
#' @param censor Detection limit for exact zeros in either Gamma model.
#'   Defaults to the family object's limit or 0.5. Zero allows strictly
#'   positive, uncensored Gamma measurements.
#' @param maxit Maximum scoring iterations per starting profile.
#' @param ncores Thread budget for sufficiently large native scan updates.
#' @return A numeric vector of rates in input scan order. Multiply by exposure
#'   for per-scan totals when the supplied spectrum sums to one.
#' @details The positive support is contiguous and the slopes of the log rate
#'   are nonincreasing in time. Weighted projection onto this cone supplies
#'   proposals for a line search against the original observation likelihood
#'   (scaled Poisson criterion for quasi-Poisson). Exact zeros in Gamma fits
#'   contribute a left-censoring probability. This is a local optimization.
#'   For multiple components [logconNMF()] alternates profile, spectrum and
#'   background updates and applies its joint rank criterion.
#' @example examples/fitLogConcaveProfile.R
#' @export
fitLogConcaveProfile <- function(y, spectrum, offset = 0, family = poisson("identity"),
                                 time = NULL, exposure = 1, weights = 1, start = NULL,
                                 dispersion = NULL, theta = NULL, censor = NULL,
                                 maxit = 4L, ncores = 1L) {
  y <- as.matrix(y)
  n <- nrow(y)
  m <- ncol(y)
  if (!is.numeric(y) || !n || !m || any(!is.finite(y)))
    stop("Supply a nonempty finite numeric observation matrix.")
  if (!is.numeric(spectrum) || length(spectrum) != m || any(!is.finite(spectrum)) ||
      any(spectrum < 0) || !is.finite(sum(spectrum)) || sum(spectrum) <= 0)
    stop("Supply one nonnegative loading per channel with positive finite total.")
  spec <- slcn_family(family, dispersion, theta)
  if (spec$kind != 0L && any(y < 0)) stop("Negative measurements require Gaussian noise.")
  if (spec$kind %in% c(1L, 3L, 5L) && any(abs(y-round(y)) > 1e-7))
    stop("Poisson and negative binomial observations must be counts.")
  if (is.null(spec$dispersion)) spec$dispersion <- 1
  if (spec$kind == 3L && (is.null(spec$theta) || spec$theta_update))
    stop("Supply a fixed positive theta for the NB2 profile.")
  if (is.null(spec$theta)) spec$theta <- 1
  if (is.null(censor)) censor <- spec$family$scalar$censor
  if (is.null(censor)) censor <- 0.5
  if (length(censor) != 1L || !is.numeric(censor) || !is.finite(censor) || censor < 0)
    stop("censor must be nonnegative.")
  if (spec$kind %in% c(4L,6L) && (any(y > 0 & y < censor) || (censor == 0 && any(y == 0))))
    stop("Encode Gamma values below the detection limit as zero and use a positive limit for zeros.")
  if (is.null(time)) time <- seq_len(n)
  if (!is.numeric(time) || length(time) != n || any(!is.finite(time)) || any(diff(time) <= 0))
    stop("time must have one strictly increasing finite coordinate per row.")
  if (!is.numeric(exposure) || !length(exposure) %in% c(1L,n) ||
      any(!is.finite(exposure)) || any(exposure <= 0)) stop("Invalid scan exposure.")
  if (is.null(start)) start <- numeric(n)
  if (!is.numeric(start) || length(start) != n || any(!is.finite(start)) || any(start < 0))
    stop("start must be a finite nonnegative profile in scan order.")
  observation <- function(x) {
    if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0))
      stop("Offsets and weights must be finite and nonnegative.")
    sinmf_observation(x,n,m)
  }
  if (length(maxit) != 1L || !is.numeric(maxit) || !is.finite(maxit) ||
      maxit < 0 || maxit != floor(maxit) || maxit > .Machine$integer.max)
    stop("maxit must be a nonnegative integer.")
  slcn_fit_profile(y, seq_len(m), 0L, observation(offset), observation(weights),
                   start, spectrum, time, rep_len(exposure,n), spec$kind,
                   spec$dispersion, spec$theta, censor, as.integer(maxit), sinmf_core_count(ncores))
}
