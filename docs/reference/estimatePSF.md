# Estimate constant or time-varying point spread functions

Estimate a log-concave PSF from a calibration standard or an overlapping
sample, and construct its sparse scan-by-scan convolution design.
Numeric retention-time names and the bundled data's `alkanes_TIC`
attributes are recognised automatically.

## Usage

``` r
estimatePSF(
  y,
  family = quasipoisson("identity"),
  time = NULL,
  standard = attr(y, "alkanes_TIC"),
  method = c("auto", "standard", "blind"),
  nwindows = 1L,
  precluster = 1L,
  dispersion = NULL,
  theta = NULL,
  weights = NULL,
  offset = NULL,
  ncores = NULL,
  control = list(),
  verbose = FALSE
)
```

## Arguments

- y:

  Numeric signal vector, scan-by-channel matrix, or a list with `Y` and
  `time`. Sparse matrices are accepted. At least nine scans are needed.

- family:

  Identity-link family supported by
  [`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md),
  including Gaussian, Poisson, quasipoisson, negative binomial,
  [`nb1()`](https://tomwenseleers.github.io/L0glm/reference/signal_families.html),
  ordinary
  [`gammaCensored()`](https://tomwenseleers.github.io/L0glm/reference/signal_families.html)
  and common-scale
  [`gammaScale()`](https://tomwenseleers.github.io/L0glm/reference/signal_families.html).
  Character names are also accepted. Gaussian observations may be
  negative.

- time:

  Increasing sample retention times. Defaults to numeric names or row
  names, or scan indices when names are absent. The standard and sample
  must use the same time units.

- standard:

  Optional nonoverlapping calibration signal, with its own
  retention-time names or `list(Y=..., time=...)`. Defaults to
  `attr(y, "alkanes_TIC")`, supplied with
  [gcmsData1](https://tomwenseleers.github.io/logconNMF/reference/gcmsData1.md)
  and
  [gcmsData2](https://tomwenseleers.github.io/logconNMF/reference/gcmsData2.md).
  A result from `estimatePSF(method="standard")` reuses its fitted
  calibration directly, constructing a sparse design on the requested
  sample grid.

- method:

  `"auto"` calibrates from an available standard and then refines on the
  sample; without a standard it uses blind sample estimation.
  `"standard"` uses calibration alone, treating `y` as the standard when
  `standard=NULL`. `"blind"` estimates directly from the overlapping
  sample. A supplied `control$kernel_est` takes precedence over
  automatic calibration.

- nwindows:

  Number of PSF anchors across the sample. `1L` gives a constant PSF.
  Multiple anchors are smoothly interpolated on the log-profile scale.

- precluster:

  `FALSE` retains the original channels, `1L` sums all channels, and an
  integer requests that many outcome groups at most. `TRUE` uses at most
  32 groups. The default uses the total-ion signal. Poisson distances
  use channel totals as size factors and native OpenMP Ward or Bregman
  clustering; see
  [`glmCluster()`](https://tomwenseleers.github.io/L0glm/reference/glmCluster.html).

- dispersion, theta:

  Shared noise parameters, as in
  [`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md).
  Omitted parameters are estimated from temporal differences and held
  fixed during PSF estimation. For NB2, `theta` is the shape parameter.

- weights, offset:

  Nonnegative scalar, scan vector or matrix matching the signal used for
  estimation. In `method="standard"` these refer to the standard; in
  `"auto"` and `"blind"` they refer to the sample. An identity-link
  offset is an additive known signal. Grouped offsets are summed.

- ncores:

  Number of native threads; defaults to the physical core count.

- control:

  Named list. `L` fixes an odd kernel length; otherwise standard support
  is inferred from its fitted peak masses and blind support defaults to
  41 scans. `kernel_est` supplies a vector or one kernel row per window.
  `maxiter=8L` and `tol=1e-3` control blind iterations; `maxiter=0L`
  returns the supplied kernel design directly. `fit` is a named list of
  additional
  [`l0glmBlindDeco()`](https://tomwenseleers.github.io/logconNMF/reference/l0glmBlindDeco.md)
  fitting options, such as `IC`. Standard controls include
  `n_peaks=NULL`, `baseline_df=NULL`, `baseline_tau=0.05`, `peak_snr=5`,
  `peak_height=NULL`, `peak_min_distance=NULL`,
  `support_fraction=0.005`, `psf_mass=0.999` and
  `support_quantile=0.95`. `linkage="ward"` controls channel clustering.
  `censor=0.5` is the Gamma detection limit, taken from a supplied
  scalar family constructor unless explicitly overridden.

- verbose:

  Print the estimation mode, aggregation and fitted PSF width.

## Value

A list containing the sparse design `X` (`X_spikes`), kernel vector or
matrix `g_hat`, `psf_kernels`, relative-time `psf_offsets`, anchor
retention times `psf_timepoints`, sample `time`, fitted `width_scans`,
`method`, noise `spec`, Gamma detection limit `censor`, aggregation
diagnostics, optional `calibration`, blind `windows`, `iterations`,
`converged` and elapsed `seconds`. The result can be passed directly as
`X=psf$X` to
[`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
or
[`l0glmFit()`](https://tomwenseleers.github.io/L0glm/reference/l0glmFit.html).
No observed matrix is retained.

## Details

**Reuse a standard calibration.** Fit a solvent-trimmed ladder TIC once
with `method="standard"`, then pass that result as `standard` for
additional samples. Its local peak kernels and their retention times are
reused; the sparse design is rebuilt for each sample's own grid. The
calibration retains the noise model used to estimate it.

**Standard calibration.** A lower-quantile spline estimates background.
Isolated calibration peaks are retained using signal-to-noise and
profile geometry. Each peak is fitted on its original observation scale:
\$\$\mu\_{ij}=b\_{ij}+o\_{ij}+h_i s_j,\qquad h_i\>0,\quad \Delta^2\log
h_i\leq0.\$\$ Native likelihood-based profile updates minimize
\\-\sum\_{ij}w\_{ij}\ell(y\_{ij};\mu\_{ij})\\ while imposing hard log
concavity. The logarithm parameterizes the positive profile; likelihood
or quasi-likelihood is evaluated on the original observation scale.
Unit-area kernels are aligned by their centroid and resampled to the
sample scan interval. Isolated-peak fitting, automatic support and
interpolation share helpers with
[`estimatePSFFromStd()`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSFFromStd.md).

**Blind refinement.**
[`l0glmBlindDeco()`](https://tomwenseleers.github.io/logconNMF/reference/l0glmBlindDeco.md)
alternates sparse nonnegative spike fitting with PSF estimation.
Completed candidates are compared using the requested information
criterion and the best checkpoint is retained. The unpenalized
background and any known offset enter both updates. With
fixed-dispersion quasipoisson counts, the Poisson likelihood weighted by
\\1/\phi\\ has the same mean-dependent objective. This permits the
native Poisson backend while preserving the quasipoisson noise scale.
For multiple windows, local kernels are estimated independently, and
\\\log g_t=(1-a)\log g_l+a\log g_r\\, with smoothstep weights
\\a=u^2(3-2u)\\, joins adjacent anchors. The returned design has
nonnegative, unit-sum, log-concave columns, including truncated boundary
columns. Kernel vectors are indexed by scans. Standard profiles are
resampled using the median sample scan interval; a constant kernel gives
shifted copies on that scan grid, independently of rounding in the
retention-time labels.

**Channel aggregation.** Poisson sums retain their model.
Shared-dispersion quasipoisson and NB1 sums retain their mean-variance
relationship, Gaussian variances add, and common-scale Gamma shapes add.
NB2 and ordinary Gamma use variance matching with observed channel
fractions \\r_j=\sum_i y\_{ij}/\sum\_{i,j\in G}y\_{ij}\\:
\$\$\theta_G=\theta/\sum\_{j\in G}r_j^2,\qquad \phi_G=\phi\sum\_{j\in
G}r_j^2.\$\$ These are generally moment-matched sum distributions. Equal
component means give exact NB2 or ordinary Gamma closure under shared
per-channel parameters. Varying channel proportions over time also
affect the variance approximation. Censored Gamma aggregation is
approximate and flagged separately. Different observation weights within
a group are averaged and flagged in the returned diagnostics; retain the
original channels for exact weighted fitting.

Calibration and sample peaks should represent the same chromatographic
conditions. Broad overloaded standard peaks can require sample
refinement; `method="standard"` keeps the calibration unchanged.
Supplied kernels are useful when a particular experimental PSF is
already known.

## References

de Rooi, J. J. and Eilers, P. H. C. (2011). Deconvolution of pulse
trains with the L0 penalty. *Analytica Chimica Acta* 705, 218–226.
[doi:10.1016/j.aca.2011.05.030](https://doi.org/10.1016/j.aca.2011.05.030)
.

## See also

[gcmsDeco()](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.html),
[`l0glmBlindDeco()`](https://tomwenseleers.github.io/logconNMF/reference/l0glmBlindDeco.md),
[`estimatePSFFromStd()`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSFFromStd.md),
[`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md),
[`glmCluster()`](https://tomwenseleers.github.io/L0glm/reference/glmCluster.html),
[`convKernel()`](https://tomwenseleers.github.io/logconNMF/reference/convKernel.md),
[`convMatrix()`](https://tomwenseleers.github.io/logconNMF/reference/convMatrix.md)

## Examples

``` r
library(L0glm)
set.seed(12)
n <- 500L
time <- seq_len(n) / 60
L <- 51L
kernel <- convKernel("bigaussian", L,
  list(sigma1 = 1.25, sigma2 = 5), perc = "mean")
X <- convMatrix(kernel, n)
A <- matrix(0, n, 3L)
A[c(60, 180, 320, 440), ] <- 30000
Y <- matrix(rpois(n * 3L, as.matrix(X %*% A) + 10), n, 3L)
rownames(Y) <- time

# Calibrate from isolated peaks with the same asymmetric instrumental response.
standard_spikes <- numeric(n)
standard_spikes[c(80, 250, 420)] <- c(90000, 120000, 100000)
standard <- list(Y = rpois(n, as.numeric(X %*% standard_spikes) + 20), time = time)
calibration <- estimatePSF(Y, time = time, family = poisson("identity"),
  standard = standard, method = "standard", nwindows = 1L, ncores = 2L,
  control = list(L = L, baseline_df = 4L, n_peaks = 3L))
dim(calibration$X)

# Blind estimation starts from a symmetric Gaussian, using summed channels.
initial <- convKernel("gaussian", L, list(sigma = 4))
blind <- estimatePSF(Y, time = time, family = poisson("identity"),
  standard = NULL, method = "blind", precluster = 1L, offset = 10,
  ncores = 2L, control = list(kernel_est = initial, maxiter = 6L))
blind$width_scans

lag <- (seq_len(L) - (L + 1L) / 2) / 60
matplot(lag, cbind(kernel, initial, calibration$g_hat, blind$g_hat),
  type = "l", lty = c(1, 2, 1, 1), lwd = 2,
  col = c("black", "grey60", "#008b9a", "#cc3344"),
  xlab = "Relative retention time (min)", ylab = "Unit-area PSF")
legend("topright", c("True bigaussian", "Initial Gaussian", "Calibrated", "Blind"),
  lty = c(1, 2, 1, 1), lwd = 2,
  col = c("black", "grey60", "#008b9a", "#cc3344"), bty = "n")

if (FALSE) { # \dontrun{
data(gcmsData2)
Y <- gcmsData2
# The alkane TIC and its time axis are taken from attr(Y, "alkanes_TIC").
psf <- estimatePSF(Y, nwindows = 3L)
B <- splines::bs(as.numeric(rownames(Y)), df = 4L, intercept = TRUE)
fit <- splogconNMF(psf$X, Y, background = B,
  family = quasipoisson("identity"), precluster = 32L)

# Calibration alone, or a single common PSF, use the same interface.
calibrated <- estimatePSF(Y, method = "standard", nwindows = 3L)
constant <- estimatePSF(Y, nwindows = 1L)
} # }
```
