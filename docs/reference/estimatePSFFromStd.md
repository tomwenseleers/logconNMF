# Estimate a time-varying PSF matrix from an alkane ladder TIC

Estimates the time-varying chromatographic point-spread function (PSF)
from the total ion chromatogram of an alkane ladder standard and returns
a sparse convolution design matrix. The baseline is estimated as a
low-degree B-spline quantile regression underapproximation. Major ladder
peaks are detected in the baseline-corrected TIC, each local peak
profile is smoothed as a non-negative, centroid-centred, log-concave
unit-sum PSF, and the local log-PSFs are interpolated or extrapolated
over retention time. The default support is inferred from the complete
local peak windows and the output scan spacing, allowing broad and
narrow peaks to share a suitable grid without imposing a fixed number of
support points.

The returned sparse matrix `X` has one column per requested output
retention time. Column `j` contains the PSF centred at
`output_timepoints[j]` and is normalised to sum to one. It can be passed
as the `X` argument to
[`l0glmSparseDeco`](https://tomwenseleers.github.io/logconNMF/reference/l0glmSparseDeco.md)
for deconvolution with a time-varying PSF.

## Usage

``` r
estimatePSFFromStd(
  y,
  family = stats::poisson(link = "identity"),
  baseline_df = 10L,
  output_timepoints = NULL,
  time = NULL,
  baseline_tau = 0.05,
  n_peaks = NULL,
  peak_min_distance = NULL,
  peak_height = NULL,
  peak_snr = 5,
  support_fraction = 0.005,
  psf_mass = 0.999,
  support_quantile = 0.95,
  kernel_half_width = NULL,
  psf_length = NULL,
  psf_nknots = NULL,
  psf_tau = 0.75,
  psf_weight_power = 0.5,
  log_floor = 1e-12,
  drop_tol = 1e-10,
  return_details = FALSE,
  verbose = FALSE
)
```

## Arguments

- y:

  Numeric vector containing the alkane ladder total ion chromatogram. If
  `time` is omitted and `names(y)` are numeric, the names are used as
  retention times.

- family:

  GLM family object with an identity link, or one of `"poisson"`,
  `"gaussian"` or `"gamma"`. The family is used to validate the intended
  non-negative identity-link noise model and to set family-aware weights
  for log-concave local PSF smoothing.

- baseline_df:

  Integer; degrees of freedom / number of spline knots for the
  low-dimensional B-spline baseline underapproximation.

- output_timepoints:

  Numeric vector of retention times for the rows and columns of the
  returned PSF matrix. Defaults to the input `time`.

- time:

  Optional numeric vector of retention times for `y`.

- baseline_tau:

  Quantile used for the B-spline baseline underapproximation. Small
  values such as `0.02` to `0.10` are useful for chromatograms with
  positive peaks.

- n_peaks:

  Optional integer number of major ladder peaks to keep. For a C9–C40
  ladder this is typically `32`. If omitted, all sufficiently prominent
  local maxima after non-maximum suppression are used.

- peak_min_distance:

  Optional minimum distance, in retention-time units, between selected
  peaks. If omitted, a conservative value is inferred from the
  chromatogram range.

- peak_height:

  Optional minimum baseline-corrected peak height. If omitted, a robust
  threshold is inferred from the baseline residuals.

- peak_snr:

  Signal-to-noise multiplier used when `peak_height` is inferred
  automatically.

- support_fraction:

  Fraction of the peak apex height used to trim local peak windows. The
  trimming is also bounded by neighbouring selected peaks.

- psf_mass:

  Fraction of each local peak's baseline-corrected signal mass used to
  infer the automatic PSF support. The default uses the central 99.9
  percent of each detected peak and avoids very long low-intensity tails
  dominating `kernel_half_width`.

- support_quantile:

  Quantile of detected local peak half-widths used to infer
  `kernel_half_width` when it is not supplied. Larger values keep more
  tail support but make the returned sparse matrix less sparse.

- kernel_half_width:

  Optional PSF half-width in retention-time units. If omitted, it is
  estimated from the detected peak windows.

- psf_length:

  Optional odd integer number of support points used for the local PSF
  grid. If omitted, it is inferred from `kernel_half_width` and the
  median spacing of `output_timepoints`.

- psf_nknots:

  Optional number of knots for the log-concave local PSF smoothing. If
  omitted, a small value is inferred from `psf_length`.

- psf_tau:

  Quantile used by the constrained log-spline smoother for each local
  PSF. Values above `0.5` fit a more upper-envelope-like log-concave
  curve and tend to preserve heavier tails.

- psf_weight_power:

  Exponent applied to the local normalised peak profile when
  constructing log-spline weights. Poisson-like weighting corresponds to
  `1`; the default `0.5` downweights tails less and usually gives more
  realistic early alkane peak tails. Use `0` for equal weights.

- log_floor:

  Positive relative floor used before interpolating local PSFs on the
  log scale.

- drop_tol:

  Relative tolerance for dropping tiny PSF entries from the sparse
  matrix.

- return_details:

  Logical; if `TRUE`, return a list with `X` and diagnostic details. If
  `FALSE`, return `X` directly with the same diagnostics attached as
  attributes.

- verbose:

  Logical; print a short summary of the fitted baseline, peak detection,
  and PSF support.

## Value

By default, a sparse `dgCMatrix` with dimensions
`length(output_timepoints) x length(output_timepoints)`. With
`return_details = TRUE`, a list with components `X`, `baseline`,
`residual`, `peak_indices`, `peak_timepoints`, `peak_heights`,
`psf_offsets`, `psf_kernels`, and `family`.

## See also

[`estimatePSF`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSF.md)
for unified standard and blind estimation,
[`l0glmSparseDeco`](https://tomwenseleers.github.io/logconNMF/reference/l0glmSparseDeco.md),
[`logConcSmooth`](https://tomwenseleers.github.io/logconNMF/reference/logConcSmooth.md),
[`convMatrix`](https://tomwenseleers.github.io/logconNMF/reference/convMatrix.md)

## Examples

``` r
library(L0glm)
data(gcmsData1)
# total ion chromatogram of C9 to C40 alkane ladder
tic <- attr(gcmsData1, "alkanes_TIC")

plot(x=as.numeric(names(attr(gcmsData1, "alkanes_TIC"))),  
           y=attr(gcmsData1, "alkanes_TIC"), 
           xlab="Retention time (min)", 
           ylab="Total ion intensity", 
           main="C9 to C40 alkane ladder standard", 
           col = "blue",
           type="l")

# Estimate the time-varying PSF with automatic support and scan-grid length.
fit <- estimatePSFFromStd(
  y = tic,
  family = poisson(identity),
  baseline_df = 10,
  output_timepoints = as.numeric(names(tic)),
  n_peaks = 32,
  return_details = TRUE
)

# fitted peaks:
matplot(fit$psf_offsets, t(fit$psf_kernels[1:32, ]),
        type = "l", lty = 1,
        xlab = "Offset from centroid (min)",
        ylab = "PSF density",
        main = "Fitted peak shapes")

# estimated time-varying PSF:
matplot(fit$X[seq(1, nrow(fit$X), by = 1), seq(1, ncol(fit$X), by = 300)], 
        type = "l", col = "blue", lty = 1,
        xlab = "Time index", ylab = "Intensity",
        main = "Estimated time-varying PSF")

# deconvolute GC/MS chromatogram using estimated time-varying PSF:
X <- estimatePSFFromStd(
  y = tic,
  family = poisson(identity),
  baseline_df = 10,
  output_timepoints = as.numeric(rownames(gcmsData1)),
  n_peaks = 32,
  return_details = FALSE
)
matplot(X, 
        type = "l", col = "blue", lty = 1,
        xlab = "Time index", ylab = "Intensity",
        main = "Estimated time-varying PSF\nfor chromatogram subsection")
system.time(fit.mc <- l0nnPois(
  X,
  gcmsData1,
  lambda = "gic",
  grouped = TRUE,
  parallelize = TRUE,
  ncores = round(parallel::detectCores(logical = FALSE))
))
# TODO: replace with l0glmFit with quasipoisson or also allow l0nnPois to estimate the dispersion

plotMat(gcmsData1, 
        xlab = "Elution time (min)",
        ylab = "Mass (m/z)",
        main = "Section of Cataglyphis iberica\nant queen cuticular GC/MS profile")
sum(Matrix::rowSums(fit.mc$coefficients)!=0)
rows = Matrix::rowSums(fit.mc$coefficients)[which(Matrix::rowSums(fit.mc$coefficients)!=0)]
nrows = rows/max(rows)
abline(v = which(Matrix::rowSums(fit.mc$coefficients)!=0),
       col = sapply(nrows, function (alpha) adjustcolor("cyan", alpha = alpha)))


hist(Matrix::rowSums(fit.mc$coefficients)[which(Matrix::rowSums(fit.mc$coefficients)!=0)])

for (i in 1:length(which(Matrix::rowSums(fit.mc$coefficients)!=0)) ) {
  plot(x = as.numeric(colnames(gcmsData1)),
       y = fit.mc$coefficients[which(Matrix::rowSums(fit.mc$coefficients)!=0)[i],],
       main = paste0("Spectrum deconvoluted compound ", as.character(i),
                     "\nRT ", rownames(gcmsData1)[which(Matrix::rowSums(fit.mc$coefficients)!=0)[1]], " min"),
       ylab = "Intensity", xlab = "Mass",
       col = "blue", type = "h")
}
```
