# Estimate the PSF given a spike train estimate via L0-penalised GLMs

Solves the *transpose* problem: given an estimated latent single- or
multichannel spike train signal \\a\\ and the vector (or matrix) of
observations \\y\\, find the PSF \\g\\ that minimises the deviance in
\\y \approx X(g)\\a\\, where \\X(g)\\ is a Toeplitz design matrix that
contains time-shifted copies of point-spread function \\g\\ (PSF, a.k.a.
convolution kernel) and errors \\\varepsilon\\ are assumed to be drawn
from an exponential family error distribution (e.g. Poisson, Gaussian,
Gamma, Inverse Gaussian). Regularisation is again done via L0glm. An
optional prior guess for the PSF can be given, in which case it is used
as the shrinkage centre in L0glm. Optional post-processing applies a
log-concave spline smooth or custom smoothing function and/or recentres
the PSF. The fitted kernel retains its positive tails; the supplied
offset represents additive background during estimation.

## Usage

``` r
l0glmDecoEstKernel(
  y,
  y_deconv,
  offset = NULL,
  kernel_est = NULL,
  L = ifelse(!is.null(kernel_est), length(kernel_est), 51),
  family = poisson(identity),
  kernel_smoothfun = function(y) logConcSmooth(y),
  perc = "mean",
  IC = "gic",
  weights = NULL,
  dispersion = NULL
)
```

## Arguments

- y:

  Numeric vector or matrix; the observed (blurred) signal. A matrix is
  treated as \\n \times m\\ outcome channels.

- y_deconv:

  Matrix or vector with the (already estimated) latent spikes.

- offset:

  Matrix or vector with a model offset (same shape as `y`), e.g. a
  previously fitted baseline/background.

- kernel_est:

  Optional prior guess for the PSF. It is L1-normalised internally and
  used as the shrinkage centre in L0glm.

- L:

  Length of the PSF to estimate (ignored when `kernel_est` is supplied).

- family:

  A GLM [family](https://rdrr.io/r/stats/family.html) object with an
  *identity* link specifying the desired error distribution (e.g.
  \`poisson(identity)\`, \`gaussian(identity)\` or \`Gamma(identity)\`).

- kernel_smoothfun:

  Optional function applied to the raw estimate (default:
  [`logConcSmooth`](https://tomwenseleers.github.io/logconNMF/reference/logConcSmooth.md)).

- perc:

  Passed to
  [`shiftKernel`](https://tomwenseleers.github.io/logconNMF/reference/convKernel.md) -
  use \`"mean"\` to recenter the centroid, or a value in \\\[0,1\]\\ for
  a percentile.

- IC:

  Information criterion for the penalty preset; see
  [information_criteria](https://tomwenseleers.github.io/L0glm/reference/information_criteria.html).

- weights:

  Optional nonnegative scalar, scan vector or matrix matching `y`.
  Applied to the kernel likelihood as well as the spike fit.

- dispersion:

  Optional fixed Gaussian or quasi-family dispersion, shared across
  channels or supplied per channel.

## Value

A list with

- `kernel`:

  The final (average) PSF estimate.

- `kernels`:

  Channel-wise estimates before averaging.

- `informative`:

  Logical vector identifying channels with nonzero spikes. Only these
  channels contribute to the average. If all spike trains are zero, the
  supplied PSF is retained, or a uniform PSF is returned.

## Details

The algorithm is inspired by the L0-penalised blind deconvolution
approach of de Rooi & Eilers (2011) for chromatographic pulse trains,
but it leverages the general machinery of L0glm. This results in the
following improvements:

- No manual tuning of the regularisation parameter is required; instead
  it can be preset to approximate information criteria such as GIC, AIC,
  BIC, EBIC or MBIC

- Both single channel or multichannel data are allowed for

- We allow the estimated convolution kernel \\g\\ to be smoothed by a
  custom function (e.g. log-concave splines)

- Any exponential-family noise model with an *identity* link can be used
  (Poisson, Gaussian, Gamma, ...)

- Shrinkage towards an optional prior guess for the PSF is allowed for,
  resulting in enhanced stability and improved convergence behaviour

- The underlying sparse design is handled by Eigen's efficient sparse
  Cholesky solver, enabling large problems to be tackled quickly

## References

de Rooi, J. J. & Eilers, P. H. C. (2011). *Deconvolution of pulse trains
with the L0 penalty*. **Analytica Chimica Acta**, *705*(1-2), 218-226.
DOI: 10.1016/j.aca.2011.05.030

## See also

[`l0glmSparseDeco`](https://tomwenseleers.github.io/logconNMF/reference/l0glmSparseDeco.md),
[`shiftKernel`](https://tomwenseleers.github.io/logconNMF/reference/convKernel.md),
[`logConcSmooth`](https://tomwenseleers.github.io/logconNMF/reference/logConcSmooth.md)

## Examples

``` r
library(L0glm)
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
                      perc   = 0)   # already centred

## --- latent sparse spike-train ------------------------------------
n        <- 10000                   # timepoints
npeaks   <- round(n / 10)           # 1000, ~10 % non-zero
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
plot(y, type = "l", ylim = c(0,max(y)), main = "Simulated blurred spike train with baseline in red")
lines(y_baseline, col="red")


## ------------------------------------------------------------------------
## 2  Deconvolution with an initial PSF guess
## ------------------------------------------------------------------------
g0 <- convKernel("gaussian", L, list(sigma = 6))  # initial guess for PSF
fit <- l0glmSparseDeco(y,
                         kernel = g0,
                         baseline = baseline,
                         family = poisson(identity), 
                         IC = "gic")
a_est <- as.vector(fit$coef_spikes)               # estimated spike train
a_baseline <- fit$coef_baseline                   # estimated baseline coefs
X_baseline <- fit$X_baseline                      # design matrix for baseline 
if (!is.null(X_baseline)) { offset <- X_baseline %*% a_baseline } else { offset <- NULL }

## ------------------------------------------------------------------------
## 3  Update PSF estimate given initial spike train estimate
## ------------------------------------------------------------------------
fit <- l0glmDecoEstKernel(y,
                            y_deconv = a_est,  # estimated spike train
                            kernel_est = g0,   # for more stable & faster convergence
                            offset = offset,   # offset for baseline
                            L = L,
                            family = poisson(identity),
                            IC = "gic")
g_est <- fit$kernel

plot(g_true, col="red",
     lwd=7, type="l",
     ylab = "Density",
     main = paste("True vs. estimated bigaussian PSF",
                  "after 1 Poisson L0glm deconvolution iteration",
                  sep = "\n"))
lines(g_est, col="blue", lwd=3)
```
