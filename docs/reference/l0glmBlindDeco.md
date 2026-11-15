# Alternating sparse L0 glm powered blind deconvolution

Runs a simple block-coordinate (alternating) optimisation: repeatedly

1.  fix the PSF and estimate spikes with `l0glmSparseDeco`;

2.  fix the spikes and estimate the PSF with `l0glmDecoEstKernel`;

until `maxiter` iterations have elapsed or the relative L1 change in the
PSF drops below `tol`. This allows one to simultaneously determine the
point spread function and spike locations of an observed blurred single-
or multichannel spike train corrupted by exponential-family noise
(Poisson, Gaussian, Gamma, Inverse Gaussian, etc.). The centering
specified by `perc` sets the time origin of the estimated PSF and
spikes. Use the same convention for simulated reference kernels when
comparing spike locations or coefficients directly. Set `perc = 0` to
retain the estimated time origin without recentering.

## Usage

``` r
l0glmBlindDeco(
  y,
  family = poisson(identity),
  kernel_est = NULL,
  L = if (is.matrix(kernel_est)) ncol(kernel_est) else if (!is.null(kernel_est))
    length(kernel_est) else 41,
  baseline = NULL,
  nknots = max(4L, floor(nrow(as.matrix(y))/500)),
  kernelfunc = function(x) convKernel("gaussian", L, list(sigma = 6)),
  kernel_smoothfun = function(y) logConcSmooth(y),
  perc = "mean",
  IC = "gic",
  grouped = TRUE,
  kernel_true = NULL,
  plot = FALSE,
  maxiter = 15,
  tol = 0.001,
  nwindows = 1L,
  ...
)
```

## Arguments

- y:

  Numeric vector or matrix; the observed (blurred) signal. A matrix is
  treated as \\n \times m\\ outcome channels.

- family:

  A GLM [family](https://rdrr.io/r/stats/family.html) object with an
  *identity* link specifying the desired error distribution (e.g.
  \`poisson(identity)\`, \`gaussian(identity)\` or \`Gamma(identity)\`).

- kernel_est:

  Initial PSF (L1-normalised internally). For `nwindows > 1`, a matrix
  with one row per window can supply local starting PSFs of equal
  length.

- L:

  Desired kernel length when \`kernel\` is \`NULL\`.

- baseline:

  Character. One of `NULL` (no baseline, default), `"intercept"` (a
  single constant column), or `"bspline"` (cubic B-spline basis). The
  baseline columns are prepended to the Toeplitz design and are *always
  included* in the model. Their fitted contribution is used as an
  `offset` in the `l0glmDecoEstKernel` kernel-update step.

- nknots:

  Integer \>= 4. Number of spline knots to use when
  `baseline = "bspline"`. With `nknots < 4` and `baseline="bspline"`
  `baseline` is changed to `intercept`.

- kernelfunc:

  A function that returns a length-\\L\\ numeric vector when called as
  \`kernelfunc(1:L)\`. Defaults to a Gaussian with sigma=6.

- kernel_smoothfun:

  Optional smoother applied to the raw PSF in each iteration (default:
  [`logConcSmooth`](https://tomwenseleers.github.io/logconNMF/reference/logConcSmooth.md)).

- perc:

  Passed to
  [`shiftKernel`](https://tomwenseleers.github.io/logconNMF/reference/convKernel.md)
  when the PSF is updated. Either '"mean"' (shift the centroid, default)
  or a numeric percentile in \[0, 1\] to shift a given percentile in the
  CDF of the estimated kernel.

- IC:

  Information criterion for the penalty preset; see
  [information_criteria](https://tomwenseleers.github.io/L0glm/reference/information_criteria.html).

- grouped:

  Logical, passed unchanged to `l0glmFit`, specifying whether to perform
  grouped L0 feature selection, i.e. best subset of group selection
  (BSGS), in case there are multiple outcome channels. Defaults to TRUE.

- kernel_true:

  Optional vector of length L with convolution kernel / PSF ground truth
  for cases where simulated data is used. Can be plotted alongside the
  estimated PSF in each iteration when `plot=TRUE`.

- plot:

  Boolean specifying whether to plot the estimated PSF in each iteration
  (alongside the true PSF if `kernel_true` is given).

- maxiter:

  Integer \>= 0. Maximum alternating steps. Default 15. Zero keeps a
  supplied `kernel_est` fixed and fits its spike train and background
  without updating the PSF.

- tol:

  Relative L1 change in the PSF that triggers early stop (default 1e-3).

- nwindows:

  Number of overlapping windows used to estimate a time-varying PSF. The
  default `1L` estimates one shared PSF. Larger values fit local PSFs,
  interpolate their log densities smoothly across scans, and refit the
  original outcomes using the resulting sparse design. For multichannel
  Poisson data, local PSF calibration uses the summed channels; the
  final fit retains every original channel.

- ...:

  Execution controls passed to the spike/background
  [`l0glmFit()`](https://tomwenseleers.github.io/L0glm/reference/l0glmFit.html)
  calls, such as `ncores` and `nthreads`. Built-in scalar families use
  native fitting automatically, with `ncores` limiting OpenMP
  outcome-channel threads. Set `ncores=1L` for serial fitting. Custom
  callbacks remain serial and use this budget for BLAS unless `nthreads`
  overrides it.

## Value

A list with components

- `a_hat`:

  Final spike estimate.

- `g_hat`:

  Final unit-sum PSF, or a matrix with one PSF per window when
  `nwindows > 1`.

- `deviance_trace`:

  Deviance value after each iteration.

- `fitted_values`, `fitted_baseline`, `coef_baseline`:

  Fitted signal, background and background coefficients on the original
  data scale.

- `X`, `X_spikes`, `X_baseline`:

  Final sparse PSF design (`X` and `X_spikes`) and background basis.

- `psf_timepoints`, `windows`:

  For time-varying fits, calibration centers in scan coordinates and
  local diagnostics.

- `deviance`, `iterations`, `converged`:

  Final deviance, number of alternating steps and PSF convergence flag.

- `best_iteration`, `selection_score`, `selection_trace`:

  For a shared PSF, the selected update, its information criterion and
  the criterion at each completed candidate fit. Local diagnostics for
  time-varying PSFs are stored in `windows`.

Background coefficients are fitted with the spikes in each iteration.
The kernel update uses the fitted background as a channel-specific
offset. Original observations are used in every fit. A final
spike/background fit ensures that the returned coefficients and deviance
use the returned PSF. Time-varying fits blend local log PSFs using
nonnegative smoothstep weights, with linear log-tail continuation when
their supports differ. This preserves log-concavity when the local PSFs
are log-concave. Edge columns are truncated using the same convention as
[`convMatrix()`](https://tomwenseleers.github.io/logconNMF/reference/convMatrix.md).
The best completed spike/background fit according to `IC` is retained
across PSF updates. `best_iteration` identifies its PSF update; zero
denotes the supplied starting kernel. `selection_trace` records the
criterion for every completed candidate fit. Each blind update compares
the full positive-tail PSF with a concentrated proposal obtained by
removing its raw median pedestal before smoothing. Both proposals are
refitted to the original observations; `IC` chooses between them. This
lets the blind search escape broad starting shapes without imposing tail
removal on a calibrated or directly estimated PSF.

## Details

The algorithm is inspired by the L0-penalised blind deconvolution
approach of de Rooi & Eilers (2011) for chromatographic pulse trains,
but it leverages the general machinery of L0glm. This results in the
following improvements:

- No manual tuning of the regularisation parameter is required; instead
  it can be preset to approximate information criteria such as GIC, AIC,
  BIC, EBIC or MBIC

- Both single channel data with L0 best subset selection or multichannel
  data with grouped L0 best subset of group selection are allowed for

- We allow the estimated convolution kernel \\g\\ to be smoothed by a
  custom function (e.g. log-concave splines)

- Any exponential-family noise model with an *identity* link can be used
  (Poisson, Gaussian, Gamma, ...)

- Estimates of the PSF are shrunk towards an optional prior guess or the
  estimate from the previous iteration and warm starts are used, which
  alongside PSF centroid centering results in improved stability and
  better convergence behaviour

- The underlying sparse designs are handled by Eigen's efficient sparse
  Cholesky solver, enabling large problems to be tackled quickly

In the non-negative-matrix-factorisation (NMF) literature this model
would be referred to as *single-atom convolutional NMF* or *sparse
convolutional dictionary learning*, whereby a one-dimensional signal (or
multichannel matrix) \\y\\ is approximated as \\W(g)\\H\\, where
\\W(g)\\ is a banded Toeplitz matrix that holds time-shifted copies of
an estimated optimal dictionary atom/PSF \\g\\, and \\H\\ is a
(potentially multi-column) matrix of sparse, non-negative activations
\\a\\. However, relative to standard single-atom convolutional NMF we
offer several extensions:

- *Optimal sparsity*: we do not fix the sparsity level of the
  activations \\a\\ in advance; instead we let L0glm perform best-subset
  or best subset of group selection through an L\\\_{0}\\ or grouped
  L\\\_{0}\\ penalty,

- *Smoothness and shape-constraints on \\g\\*: we allow the estimated
  convolution kernel atom \\g\\ to be smoothed by a custom function
  (e.g. log-concave splines)

- *General noise models*: any exponential-family likelihood with an
  *identity* link is supported (Poisson, Gaussian, Gamma,
  inverse-Gaussian, ...), instead of only Gaussian (Euclidean-distance
  NMF) or Poisson (NMF with KL divergence)

- *Speed*: the Toeplitz designs in the alternating updates are solved by
  sparse-Cholesky linear algebra rather than multiplicative updates,
  giving orders-of-magnitude speed-ups on long signals

Methodologically this could be referred to as *generalised / deviance
matrix factorisation* (Kidzinski et al. 2022; Castiglione et al. 2024;
Wang & Carvalho 2023), specialised to the case where \\W\\ is a Toeplitz
single-atom convolutional dictionary matrix and where the matrix
factorisation is solved via alternating block-coordinate descent, with
each update being powered by a L\\\_{0}\\ or grouped
L\\\_{0}\\-penalised generalised linear model (GLM) fit via adaptive
ridge regularised IRLS.

A classical baseline for Poisson deconvolution is the Richardson-Lucy
algorithm (Richardson 1972; Lucy 1974). In modern terms, the blind
deconvolution version of Richardson-Lucy equals Kullback-Leibler NMF
with multiplicative updates (Lee & Seung 1999) and factor matrix \\W\\
constrained to be a convolution matrix. Empirically, `l0glmBlindDeco()`
has much better support recovery, attains much lower reconstruction
error and is dramatically faster thanks to (i) closed-form sparse linear
solves, (ii) the automatic selection of the penalty parameter through
information-criterion presets, (iii) semi-parametric shape constraints
on the estimated convolution kernel and (iv) the possibility to exploit
grouped-sparsity when several channels are deconvolved jointly.

## References

- Richardson, W. H. (1972)\
  Bayesian-based iterative method of image restoration.\
  *J. Opt. Soc. Am.* 62 (1): 55-59.

- Lucy, L. B. (1974)\
  An iterative technique for the rectification of observed
  distributions.\
  *Astron. J.* 79: 745-754.

- Smaragdis, P. (2004)\
  Non-negative matrix factor deconvolution; extraction of multiple sound
  sources from monophonic inputs.\
  *Proceedings of ICA 2004*.

- Lee, D. D., & Seung, H. S. (1999)\
  Learning the parts of objects by non-negative matrix factorization.\
  *Nature* 401 (6755): 788-791.

- de Rooi, J., & Eilers, P. H. C. (2011)\
  Deconvolution of pulse trains with the L0 penalty.\
  *Analytica Chimica Acta* 705 (1-2): 218-226.

- Kidzinski, L., Hui, F. K. C., Warton, D. I. and Hastie, T. (2022)\
  Generalized Matrix Factorization: efficient algorithms for fitting
  generalized linear latent variable models to large data arrays.\
  *Journal of Machine Learning Research* 23(291): 1–29.
  <https://www.jmlr.org/papers/v23/20-1104.html>.

- Castiglione, C., Segers, A., Clement, L. and Risso, D. (2026)\
  Stochastic gradient descent estimation of generalized matrix
  factorization models with application to single-cell RNA sequencing
  data.\
  *Biostatistics* 27(1), kxag010.
  [doi:10.1093/biostatistics/kxag010](https://doi.org/10.1093/biostatistics/kxag010)
  .

- Wang, L. and Carvalho, L. E. (2023)\
  Deviance matrix factorization.\
  *Electronic Journal of Statistics* 17(2): 3762–3810.
  [doi:10.1214/23-EJS2174](https://doi.org/10.1214/23-EJS2174) .

## See also

[`l0glmSparseDeco`](https://tomwenseleers.github.io/logconNMF/reference/l0glmSparseDeco.md),
[`l0glmDecoEstKernel`](https://tomwenseleers.github.io/logconNMF/reference/l0glmDecoEstKernel.md)

## Examples

``` r
library(L0glm)
# \donttest{
  
## ------------------------------------------------------------------------
## 1  Simulate a blurred spike-train with a BIGAUSSIAN PSF & Poisson noise
## ------------------------------------------------------------------------
library(Matrix)
library(splines)
set.seed(123)
family <- poisson(identity)

## --- true point-spread function (PSF) -----------------------------
L       <- 51                       # kernel length
sigma1  <- 1.25                     # left-hand sigma
sigma2  <- 5                        # right-hand sigma
g_true  <- convKernel("bigaussian", L,
                      params = list(sigma1 = sigma1,
                                    sigma2 = sigma2),
                      perc   = "mean") # match the fitted PSF's centroid convention

## --- latent sparse spike-train ------------------------------------
n        <- 10000L                  # timepoints
npeaks   <- 1000L                   # nonzero peaks
idx      <- sample(n, npeaks)
amps     <- sample(1000:100000, npeaks)
a_true   <- numeric(n);  a_true[idx] <- amps   # true spikes

## --- baseline -----------------------------------------------------
baseline <- "intercept" # or "bspline" or NULL
if (!is.null(baseline)) {if (baseline=="intercept") {
  X_baseline <- Matrix::Matrix(1, n, 1, sparse = TRUE)
  a_baseline <- 5000
} else if (baseline=="bspline") {
  nknots <- max(4L, floor(n/500))
  X_baseline <- splines::bs(seq_len(n), df = nknots, degree = 3)
  X_baseline <- Matrix::Matrix(X_baseline, sparse = TRUE)
  a_baseline <- runif(nknots, min = 2000, max = 6000)
}} else { X_baseline <- NULL }

## --- generate blurred, noisy observations -------------------------
X_sim  <- convMatrix(g_true, n)
if (!is.null(X_baseline)) { X_full <- cbind(X_baseline, X_sim)
a_full <- c(a_baseline, a_true)
y_baseline <- as.vector(X_baseline %*% a_baseline) # noiseless baseline
} else {
  X_full <- X_sim
  a_full <- a_true
  y_baseline <- rep(0, n) }

y_true <- as.vector(X_full %*% a_full)    # noiseless convolution
y      <- rpois(n, y_true)                # with Poisson noise

# quick sanity plot
plot(y, type = "l", ylim = c(0,max(y)),
     main = "Simulated blurred spike train with baseline in red")
lines(y_baseline, col="red")


## ------------------------------------------------------------------------
## 2  Blind deconvolution starting from intial Gaussian PSF guess
## ------------------------------------------------------------------------
g0 <- convKernel("gaussian", L, list(sigma = 6))  # initial guess for PSF

fit <- l0glmBlindDeco(y,
                        family = poisson(identity),
                        L = L,                 # PSF kernel length
                        baseline = baseline,
                        kernel_est = g0,
                        kernel_smoothfun = function(y) {
                          logConcSmooth(y, nknots = round(L / 5))
                        },
                        perc = "mean",         # use centroid centering
                        kernel_true = g_true,  # for plots
                        IC = "gic",
                        plot = TRUE,
                        tol = 1e-3)









plot(fit$deviance_trace, type="l", log="y",
     xlab="Iteration", ylab="Deviance",
     main="Poisson L0 glm blind deconvolution\nconvergence behaviour")

a_est = as.vector(fit$a_hat)   # estimated spikes
g_est = fit$g_hat              # estimated PSF

## ---------- summary statistics ------------------------------------
n_true <- sum(a_true != 0)
n_est  <- sum(a_est  != 0)
tp     <- sum(a_true != 0 & a_est != 0)
precision <- tp / n_est
recall    <- tp / n_true
F1        <- 2 * precision * recall / (precision + recall)
R2p <- 100 * summary(lm(a_true ~ 0 + a_est))$r.squared

## ---------- scatter-plot with true vs. estimated peak areas -------
main_expr <- bquote(
  atop("Poisson L0 glm blind deconvolution",
       atop("true vs. estimated peak areas",
            italic(R)^2 == .(sprintf("%.2f", R2p)) * "%" ~ "|" ~
              "peaks:" ~ .(n_true) ~ "(true)," ~ .(n_est) ~ "(est.)" ~ "|" ~
              "F1 =" ~ .(sprintf("%.2f", 100 * F1)) * "%")))
plot(a_est, a_true,
     xlim = range(a_est, a_true),
     ylim = range(a_est, a_true),
     pch  = 16, col = "blue",
     xlab = "estimated value",
     ylab = "true value",
     main = main_expr)


# Multichannel fits use an individual background offset for each channel.
n_multi <- 500L
keep_multi <- seq_len(n_multi)
Y <- cbind(y[keep_multi],
           rpois(n_multi, 0.7 * y_true[keep_multi]),
           rpois(n_multi, 1.3 * y_true[keep_multi]))
multi <- l0glmBlindDeco(Y, family = poisson(identity), kernel_est = g0,
                         baseline = "bspline", nknots = 6L,
                         maxiter = 3L, plot = FALSE)
dim(multi$a_hat)
#> [1] 500   3
multi$deviance
#> [1] 187483.5

# Estimate smoothly varying PSFs in ten overlapping windows.
# The returned X can also be supplied to splogconNMF().
varying <- l0glmBlindDeco(y, family = poisson(identity), kernel_est = g0,
                           baseline = "bspline", nknots = 8L,
                           nwindows = 10L, maxiter = 6L)
dim(varying$g_hat)  # 10 local PSFs
#> [1] 10 51
dim(varying$X)      # sparse scan-by-scan design
#> [1] 10000 10000
# }
```
