# Fit a log-concave factorization with automatic rank selection

Estimate ordered, nonnegative log-concave profiles and channel loadings
from one signal or several related samples. The information criterion
selects the number of components. Simultaneous fits share spectra while
retaining the original scan grids and estimating a separate profile and
area for each occurrence.

## Usage

``` r
logconNMF(
  y,
  X = NULL,
  family = poisson("identity"),
  lambda = "gic",
  background = NULL,
  dispersion = NULL,
  theta = NULL,
  precluster = "auto",
  weights = NULL,
  offset = NULL,
  correlation = NULL,
  ncores = NULL,
  control = list(),
  start = NULL,
  path = NULL,
  align = TRUE,
  method = c("cosine", "hellinger"),
  verbose = FALSE
)

# S3 method for class 'logconNMF_fit'
coef(object, sample = NULL, normalize = c("none", "l1"), ...)

# S3 method for class 'logconNMF_fit'
fitted(object, sample = 1L, ...)

# S3 method for class 'simultaneousNMF'
plot(
  x,
  type = c("areas", "spectra", "profiles"),
  sample = NULL,
  compounds = seq_len(min(6L, nrow(x$P))),
  ...
)

sampleNMF(object, sample = 1L)
```

## Arguments

- y:

  A numeric scans-by-channels matrix or vector for one sample, or a list
  of at least two sample descriptions for a simultaneous fit. Each
  description contains `y` and increasing `time`, with optional `X`,
  `background`, `fit`, `weights`, `offset`, `exposure`, `dispersion`,
  `theta`, `sample_id`, `group`, `small_width` and `edges`. Supply
  observations, weights, offsets or initial fits as RDS paths to read
  them from disk. Samples share channel names and order. List names
  supply sample IDs when the descriptions omit `sample_id`.

- X:

  Candidate profile matrix for a single sample, such as the sparse
  design from
  [`estimatePSF()`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSF.md).
  In a simultaneous fit, supply it in each sample description or let the
  simultaneous initializer estimate it.

- family:

  Identity-link noise family: `gaussian`, `poisson`, `quasipoisson`,
  [`mgcv::nb`](https://rdrr.io/pkg/mgcv/man/negbin.html), `Gamma`,
  [`nb1()`](https://tomwenseleers.github.io/L0glm/reference/signal_families.html),
  [`gammaScale()`](https://tomwenseleers.github.io/L0glm/reference/signal_families.html)
  or
  [`gammaCensored()`](https://tomwenseleers.github.io/L0glm/reference/signal_families.html).
  The corresponding names, including `"nb"`, `"nb1"` and
  `"Gamma_scale"`, are also accepted. Count likelihoods require integer
  counts; quasipoisson accepts nonnegative intensities. Gaussian
  observations may be negative; their fitted latent signals remain
  nonnegative.

- lambda:

  An information criterion (`"gic"`, `"aic"`, `"bic"`, `"ebic"` or
  `"mbic"`) or a nonnegative numeric penalty per selection degree of
  freedom, before any temporal information adjustment. With `m`
  channels, a single component incurs `m * lambda`; in a simultaneous
  fit a shared spectrum incurs `(m - 1) * lambda` and an occurrence
  incurs `lambda`.

- background:

  Optional nonnegative scans-by-basis matrix, for example
  `splines::bs(time, df=6, intercept=TRUE)`. Coefficients are
  nonnegative and unpenalized. `NULL` fits without a background.

- dispersion:

  Known common Gaussian variance, quasipoisson scale, Gamma dispersion,
  NB1 variance-to-mean ratio or common Gamma scale. `NULL` estimates a
  pooled value from temporal differences and freezes it. Parameters
  supplied in a
  [signal_families](https://tomwenseleers.github.io/L0glm/reference/signal_families.html)
  object are used when this argument is omitted. Poisson and NB2 use
  unit dispersion.

- theta:

  NB2 shape, shared by all original channels. `NULL` estimates and
  freezes a temporal moment estimate; a positive value fixes it; a
  negative value initializes at its absolute value and enables pooled
  likelihood updates. A fixed shape in the supplied NB family is
  respected.

- precluster:

  Maximum channel groups or a saved
  [`channelGroups()`](https://tomwenseleers.github.io/logconNMF/reference/channelGroups.md)
  map. `"auto"` uses the original channels for a single sample and up to
  32 common groups for a simultaneous fit. `FALSE` or `NULL` disables
  aggregation. Final fits use every original channel.

- weights:

  Nonnegative likelihood weights: a scalar, row vector or matrix
  matching `y`. `NULL` uses one. Zero weights omit observations.

- offset:

  Known nonnegative mean contribution, in the same supported formats as
  `weights`. `NULL` uses zero. It is included on the response scale
  because the link is identity.

- correlation:

  Lag-one working information adjustment. `NULL` uses the fitter's
  default (`"auto"` for one sample and `"none"` for simultaneous fits),
  or `control$correlation` if supplied for a simultaneous fit.

- ncores:

  Total core budget. `NULL` uses the detected physical core count.
  Simultaneous fits divide that budget between sample workers and their
  native threads, with one BLAS thread per worker.

- control:

  Named controls. See
  [`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  for one matrix. For a simultaneous fit, `cycles=10L` and
  `initial_cycles=5L` set refinement limits, `rounds=3L` sets the number
  of grouping rounds, and `q=64L` sets the local profile update
  resolution. The default `rt_tolerance=NULL` and
  `alignment_window=NULL` derive matching windows from the narrow-peak
  widths. `ridge=TRUE` checks the profile ridge across samples. The
  remaining defaults are `censor=0.5`, `correlation="none"`,
  `workers=NULL` and `cache=TRUE`. `initial=list()` passes controls to
  the independent initial fits. `workers` limits concurrent samples
  within `ncores`. The cache stores aggregated observations in memory.

- start:

  Optional list `H`, `S` and `alpha`, with original dimensions. `H` has
  scans in rows and compounds in columns; `S` has compounds in rows and
  all original channels in columns. `alpha` contains background
  coefficients and defaults to zero. Profiles are normalized with the
  inverse change applied to spectra. A restart bypasses initial
  screening and preclustering; nonnegative spectral entries can reopen
  during refitting.

- path:

  New directory for the sample store and fitted factors in a
  simultaneous fit. `NULL` uses a temporary directory.

- align:

  Estimate monotone time warps for a simultaneous fit.

- method:

  Spectral distance for matching simultaneous samples, either `"cosine"`
  or `"hellinger"`.

- verbose:

  Print stage timings, objective changes and compound counts.

- object:

  A fitted model returned by `logconNMF()`. `sampleNMF()` takes a
  simultaneous model.

- sample:

  Sample index or ID. For a single matrix use `1L` or `"sample1"`. In
  [`coef()`](https://rdrr.io/r/stats/coef.html), `NULL` returns the
  common spectra for a simultaneous fit. In
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), `NULL`
  displays all samples.

- normalize:

  `"none"` returns amplitude-scaled loadings. `"l1"` returns unit-sum
  rows.

- ...:

  Further arguments passed to the underlying extractor.

- x:

  A fitted simultaneous model for plotting.

- type:

  Plot sample areas, common spectra or aligned compound profiles.

- compounds:

  Compound indices to display in a simultaneous plot.

## Value

A fitted object with estimated `rank`, the number of single-sample
components or shared components. Use
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
in either mode. [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
returns one sample on its original grid. A simultaneous fit includes
consensus spectra `P`, sample-by-compound `areas`, `compounds`,
`alignment`, `paths`, `store`, `channel_groups` and `parallel` thread
allocation.

## Details

**Rank-selection objective.** The fit approximately minimizes
\$\$-2\sum_s \ell\_{\nu_s}(M_s;Y_s) + \kappa\lambda d.\$\$ The weighted
log-likelihood \\\ell\_{\nu_s}\\ uses the fitted mean \\M_s\\ on the
original channels of sample \\s\\. For quasipoisson it is the Poisson
log score divided by the sample's estimated or supplied dispersion. The
optional working temporal adjustment is \\\kappa\\, equal to one under
independence. A numeric `lambda` sets the penalty per selection degree
of freedom. The GIC preset uses \\\log(p)\log\log(n)\\, with the PSF
dictionary width and number of scans for a single matrix. For a sample
list the preset uses \\p=n=\sum_s n_s\\.

For one matrix with \\m\\ channels and \\K\\ components, \\d=mK\\. For a
simultaneous model, \\d=(m-1)G+L\\, where \\G\\ is the number of shared
unit-sum spectra and \\L\\ counts the detected sample–component
occurrences. Each spectrum contributes \\m-1\\ coordinates and each
occurrence contributes one amplitude. Setting \\G=L=K\\ gives the
single-sample count. Profiles have contiguous support and a concave
logarithm on the positive interval. Profile shapes and smooth
backgrounds carry no additional rank charge.

The simultaneous `objective` and reported information criteria such as
[`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
also count estimated noise parameters. The single-sample `objective`
drops that constant charge for a fixed noise specification. Local
refinement and fusion, deletion or addition of candidate components
provide the approximate rank search. Final single-sample refinement also
tries broader local starts for nearly coincident profiles, with
acceptance checked on every original channel against the current fit and
a refined control. This is enabled by default. For one matrix,
`control=list(profile_restarts=FALSE)` disables the pass; simultaneous
fits pass this option through `control$initial`.

For a simultaneous fit, put `background`, `weights`, `offset` and any
initial `fit` in the sample descriptions. `path = NULL` uses a new
temporary sample store. Supply a persistent path when saving a
simultaneous fit for later use, and retain the store with the fitted
object.

Numeric penalties in `logconNMF()` have the same scale in both modes. To
reproduce a single-sample
[`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
call using a numeric compound penalty, divide that penalty by the number
of channels.

## Simultaneous models

Ignoring backgrounds, exposures and offsets, the sample means have the
shared-factor form \\M_s = Q_s P\\, with one common spectral factor
\\P\\ and sample-specific profile/area factors \\Q_s\\. This follows the
simultaneous-NMF terminology of [Cichocki et al.
(2009)](https://doi.org/10.1002/9780470747278) and [Lee et al.
(2012)](https://doi.org/10.1371/journal.pone.0048238). It does not
introduce a separate cross-modality coupling matrix.

Samples can supply separate background bases, known exposures and noise
parameters. With a time warp \\u=w_s(t)\\, binwise Jacobian \\J\\, and
exposure \\e\\, the fitted mean on the original grid is
\$\$\mu\_{sij}=e\_{si}\\J\_{si}\sum_g a\_{sg}h\_{sg}(u\_{si})p\_{gj}
+(B_sA_s)\_{ij}\\+o\_{sij}.\$\$ Spectra sum to one and the profiles
integrate to one on the aligned grid. Unit exposure models the supplied
intensities. Use known integration durations as exposures for integrated
counts. Scan spacing and `small_width` use the same time units.

Reciprocal spectral matches provide leading-edge alignment anchors.
Compatible occurrences share spectra after checking overlap and profile
shape in area order. At most one occurrence of a compound enters each
sample. Local profile, spectral and background updates compare the full
original-channel likelihood. Undetected occurrences have `NA` areas.
Noise scales are held fixed when comparing groupings. A negative NB2
`theta` allows bounded, likelihood-improving shape updates.

The common channel map is learned from samples aggregated on common time
bins. The fitter retains only the aggregated observations in its shared
cache and reads the original channels in blocks from each sample store.
Keep the store with the saved fitted object for reconstruction and
plots.

## Extract one sample

`sampleNMF(object, sample)` returns a single-sample factorization view.
The view contains unit-sum original-grid profiles in `H`, matching
amplitude-scaled spectra in `S`, its fitted background, integrated areas
in `integrated_area`, and consensus IDs in `compound_id`. Alignment and
exposure are included in the profiles and spectra. Use
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
on the returned view. The view carries no observed matrix.

## Examples

``` r
X <- convMatrix(convKernel("gaussian", 9L), 50L)
signal <- cbind(700 * X[, 16L] + 350 * X[, 35L],
                150 * X[, 16L] + 500 * X[, 35L])
set.seed(12)
Y <- matrix(rpois(length(signal), signal), nrow(signal))
fit <- logconNMF(Y, X = X, ncores = 1L)
fit$rank
coef(fit, normalize = "l1")

second <- matrix(rpois(length(signal), 1.4 * signal), nrow(signal))
samples <- list(first = list(y = Y, time = seq_len(nrow(Y)), X = X),
  second = list(y = second, time = seq_len(nrow(Y)), X = X))
joint <- logconNMF(samples, ncores = 1L)
joint$rank
joint$areas
fitted(joint, sample = "second")[1:3, ]
second_fit <- sampleNMF(joint, sample = "second")
second_fit$compound_id
coef(second_fit)
```
