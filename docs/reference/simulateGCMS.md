# Simulate synthetic GC-MS chromatograms with ground-truth signal + background

\`simulateGCMS()\` draws \\k\\ compounds from an EI reference library
\`S\` (rows are compounds and columns are m/z values, e.g. 41–298) and
generates an \\n \times m\\ chromatogram \\D\\ whose signal is the sum
of \\k\\ Gaussian peaks, \\W\_{\mathrm{true}}H\_{\mathrm{true}}\\, and a
smooth, monotonically increasing B-spline background
\\W\_{\mathrm{background}}\\b\\H\_{\mathrm{background}}\\, where \\b\\
is the background scale.

## Usage

``` r
simulateGCMS(
  n = 10000,
  k = 500,
  S = NULL,
  peak_width_range = c(3, 5),
  mean_beta = 10000,
  sd_logbeta = 1,
  width_height_coef = 3,
  rt_height_coef = 0.02,
  sample_nr = 1,
  rt_drift_coef = 0.001,
  bg_scale = mean_beta/10,
  family = c("poisson", "gaussian"),
  sd_noise = if (match.arg(family) == "gaussian") 1 else NULL,
  seed_compounds = 123,
  seed_peaks = 456
)
```

## Arguments

- n:

  Integer \>= 1. Number of time points (rows).

- k:

  Integer \>= 1. Number of compounds/peaks.

- S:

  Spectral library. If \`NULL\` the example data
  \`data(simulatedSpectra)\` (2000 spectra x 258 m/z values) is loaded.

- peak_width_range:

  Length-2 numeric. Uniform base range for Gaussian peak width sigma
  (\*before\* any height-coupling is applied and on time point scale).

- mean_beta, sd_logbeta:

  Parameters of the log-normal height distribution.

- width_height_coef:

  Numeric \>= 0. Multiplicative widening: \$\$\sigma_j \leftarrow
  \sigma_j\bigl\[1+\texttt{coef}\\h_j\bigr\]\$\$ where
  \\h_j=(\beta_j-\min\beta)/(\max\beta-\min\beta)\\.

- rt_height_coef:

  Numeric. Positive values push high peaks to slightly larger RTs:
  \$\$c_j \leftarrow c_j
  \bigl\[1+\texttt{rt\\height\\coef}\\h_j\bigr\]\$\$

- sample_nr:

  Integer \>= 1. Injection/order number.

- rt_drift_coef:

  Numeric. Additive drift in \*time points\* per injection: \\c_j
  \leftarrow c_j + \texttt{rt\\drift\\coef}\\ (\texttt{sample\\nr}-1)\\.

- bg_scale:

  Numeric \>= 0. Multiplicative factor applied to the whole B-spline
  background.

- family:

  \`"poisson"\` (default) or \`"gaussian"\`.

- sd_noise:

  s.d. of additive Gaussian noise (\`family="gaussian"\`).

- seed_compounds:

  Integer seed controlling which library spectra are drawn.

- seed_peaks:

  Integer seed for all other random quantities (heights, widths, noise,
  ...).

## Value

A list with matrices:

- D:

  Observed chromatogram (signal + background + noise) (dense matrix).

- D_true:

  Noiseless signal
  \\W\_{\mathrm{true}}H\_{\mathrm{true}}+W\_{\mathrm{background}}H\_{\mathrm{background}}\\
  (dense matrix).

- W_true, H_true:

  Elution profiles and EI spectra of the \*k\* peaks (sparse matrices).

- W_background, H_background:

  Smooth B-spline background bases and spectra (sparse matrices).

- k_background:

  Number of B-spline basis columns (\`max(2, n/1000)\`).

- centres, widths, heights:

  Peak parameters after all applied height/RT/drift adjustments.

- widths_base:

  Base compound width before intensity-based broadening.

- compound_idx:

  Row indices of \`S\` that were selected.

- sample_nr, rt_drift_coef, width_height_coef, rt_height_coef:

  Echoed.

## Details

- A *compound-specific base width* \\\sigma\_{0}\\ is drawn once when
  the compound is first selected from the library.

- In subsequent simulated chromatograms the same compound retains the
  same \\\sigma\_{0}\\, but the *effective* width becomes \$\$\sigma
  \\=\\ \sigma\_{0}\Bigl(1 +
  \texttt{width\\height\\coef}\\h\_{\mathrm{norm}}\Bigr),\$\$ so it can
  still broaden or narrow in proportion to its (possibly resampled) peak
  height.

Gaussian peak width and RT can be coupled to peak height; an additional
\*chromatographic drift\* proportional to the injection number
\`sample_nr\` is also included. Independent random seeds allow you to
keep the same compounds while resampling peak heights (or vice-versa).

## Reference library

\`S\` must contain - column names (character) for the m/z channels (e.g.
\`"41":"298"\`), and - attribute \`attr(S, "RT")\` for the retention
times \*(numeric, unit: mins)\*.

## Author

Tom Wenseleers <tom.wenseleers@kuleuven.be> \\ Laboratory of
Socioecology and Social Evolution, KU Leuven \\
<https://bio.kuleuven.be/eeb/tw/research> \\ ORCID:
[0000-0002-1434-861X](https://orcid.org/0000-0002-1434-861X)

## Examples

``` r
library(L0glm)
data(simulatedSpectra)
sim <- simulateGCMS(n = 10000,
                     k = 500,
                     S = simulatedSpectra,
                     peak_width_range = c(3, 5),
                     mean_beta  = 1e4,
                     sd_logbeta = 1,
                     width_height_coef = 3,
                     rt_height_coef    = 0.02,
                     sample_nr = 1,
                     rt_drift_coef = 0.001,
                     bg_scale = 100,
                     family = "poisson")
plotMat(sim$D,
        xlab = "Time (min)", 
        ylab = "Mass (m/z)", 
        main = "Simulated GC/MS chromatogram")
```
