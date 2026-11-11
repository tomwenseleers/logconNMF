# Sparse log-concave nonnegative matrix factorization

Recover a parsimonious set of compounds from multichannel signals, with
nonnegative spectra, hard log-concave elution profiles and a smooth
background. Fast spectral fusion and local native updates accommodate
broad, asymmetric and shark-fin peaks while retaining all original
outcome channels.

## Usage

``` r
# S3 method for class 'splogconNMF'
residuals(
  object,
  y = object$y,
  type = c("pearson", "response"),
  rows = seq_len(object$n),
  channels = seq_len(ncol(object$S)),
  mu = NULL,
  ...
)

splogconNMF(
  X,
  y,
  family = poisson("identity"),
  background = NULL,
  lambda = "gic",
  dispersion = NULL,
  theta = NULL,
  precluster = NULL,
  weights = NULL,
  offset = NULL,
  correlation = c("auto", "none", "lag1"),
  ncores = NULL,
  control = list(),
  start = NULL,
  verbose = FALSE
)

# S3 method for class 'splogconNMF'
fitted(object, ...)

# S3 method for class 'splogconNMF'
nobs(object, ...)

# S3 method for class 'splogconNMF'
coef(object, ...)

# S3 method for class 'splogconNMF'
logLik(object, ...)

# S3 method for class 'splogconNMF'
AIC(object, ..., k = 2)

# S3 method for class 'splogconNMF'
BIC(object, ...)

# S3 method for class 'splogconNMF'
GIC(object, ...)

# S3 method for class 'splogconNMF'
EBIC(object, ...)

# S3 method for class 'splogconNMF'
MBIC(object, ...)

# S3 method for class 'splogconNMF'
print(x, ...)
```

## Arguments

- object, x:

  A fitted `splogconNMF` object.

- y:

  Numeric scans-by-channels matrix or a single-channel vector.

- type:

  Residual scale: Pearson or response.

- rows, channels:

  Optional row and channel indices for residual extraction.

- mu:

  Optional fitted means for the selected block.

- ...:

  Further method arguments.

- X:

  Nonnegative scans-by-candidates matrix, typically shifted copies of a
  point-spread function from
  [`convMatrix()`](https://tomwenseleers.github.io/logconNMF/reference/convMatrix.md)
  or
  [`l0glmBlindDeco()`](https://tomwenseleers.github.io/logconNMF/reference/l0glmBlindDeco.md).
  Columns must have contiguous log-concave positive support. Dense and
  sparse inputs are accepted and stored sparsely; columns are normalized
  to unit sum.

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

- background:

  Optional nonnegative scans-by-basis matrix, for example
  `splines::bs(time, df=6, intercept=TRUE)`. Coefficients are
  nonnegative and unpenalized. `NULL` fits without a background.

- lambda:

  `"gic"` (default), `"aic"`, `"bic"`, `"ebic"`, `"mbic"`, or a
  nonnegative numeric compound penalty. See Details for its scale.

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

  A saved
  [`channelGroups()`](https://tomwenseleers.github.io/logconNMF/reference/channelGroups.md)
  map shares the same aggregation across analyses. `NULL` or `FALSE`
  uses every channel throughout. `TRUE` uses up to 32 groups; an integer
  supplies another maximum. Size-adjusted Poisson clustering,
  aggregation, selection and fusion operate in that reduced space,
  followed by refinement on all original channels.

- weights:

  Nonnegative likelihood weights: a scalar, row vector or matrix
  matching `y`. `NULL` uses one. Zero weights omit observations.

- offset:

  Known nonnegative mean contribution, in the same supported formats as
  `weights`. `NULL` uses zero. It is included on the response scale
  because the link is identity.

- correlation:

  `"auto"` or `"lag1"` estimates an optional lag-one working adjustment
  to the selection penalty. `"none"` uses independent observations; a
  number between zero and 0.49 supplies the lag-one correlation.

- ncores:

  OpenMP workers. `NULL` uses the detected physical cores.

- control:

  Named list: `broad_repeats=1L`, `full_sweeps=1L`, `final_polishes=1L`,
  `native_sweeps=6L`, `block_sweeps=6L`, `psf_width=NULL`,
  `screen_factor=1`, `censor=0.5`, `linkage="ward"`,
  `joint_trials=NULL`, `profile_cycles=4L`, `fusion_passes=0L`,
  `profile_fusion=TRUE`, `profile_restarts=TRUE`. `joint_trials` uses
  direct Poisson fusion for Poisson data and joint local refinement for
  the other families. `psf_width` supplies the typical FWHM in scans.
  When omitted, local updates use the middle design column and
  broad-profile fusion uses the median of up to 21 representative design
  columns. `screen_factor` scales the initial selection penalty.
  `censor` is the Gamma detection limit: recorded zero intensities are
  left-censored below this value; positive values must be at least the
  limit. Zero requests an uncensored Gamma model. A supplied Gamma
  family constructor provides its own detection limit unless overridden
  here. `linkage` may be `"ward"` or `"bregman"`. `joint_trials` allows
  overlapping compounds to adjust together during local fusion
  refinement. `profile_cycles` bounds the smooth shape/background
  initialization; positive values also enable the final
  likelihood-checked local restarts. `profile_restarts=FALSE` disables
  those alternative starting profiles. positive values also enable fine
  native profile/spectrum refinement. Zero skips these shape-polishing
  cycles. `fusion_passes` optionally adds broad-anchor and
  compatible-tail searches during initialization. `profile_fusion`
  enables algebraic profile fusion, retained contained-fragment
  proposals and a final broad-fragment pass on every original channel.
  Broad fusions are compared with an equally refined unmerged local
  model. The entire support of every updated component is included in
  the local fit, with the remaining components and background held in an
  offset. `original_cleanup=TRUE` enables a final original-channel
  deletion, local fusion and residual-discovery pass. Set it to `FALSE`
  to skip that pass.

- start:

  Optional list `H`, `S` and `alpha`, with original dimensions. `H` has
  scans in rows and compounds in columns; `S` has compounds in rows and
  all original channels in columns. `alpha` contains background
  coefficients and defaults to zero. Profiles are normalized with the
  inverse change applied to spectra. A restart bypasses initial
  screening and preclustering; nonnegative spectral entries can reopen
  during refitting.

- verbose:

  Print stage timings, objective changes and compound counts.

- k:

  Complexity multiplier for
  [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html).

## Value

An object of class `splogconNMF`, containing sparse unit-sum profiles
`H`, spectra `S`, background coefficients `alpha`, the background basis,
fitted means `mu`, `rank`, resolved noise settings `spec`, original
data, offsets and weights, `objective`, `loglik`, stage ledger
`history`, `ledger_error`, `initialization`, aggregation,
`profile_fusion` and `contained_fusion` and `profile_starts`
diagnostics, controls and timings `seconds`.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns spectra and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) the full mean
matrix. [`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html),
[`EBIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
and
[`MBIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
score the returned model. Columns of `H` and the matching rows of `S`
are ordered by increasing apex scan, also returned in `apex_scan`. Their
labels identify the same compounds.

## Details

**Generalized matrix factorization.** The model is an identity-link,
nonnegative generalized matrix factorization (GMF), with a smooth
background, log-concave profiles and automatic rank selection. See
Kidzinski et al. (2022) and Castiglione et al. (2026). The function
retains the NMF name; [the mathematical
article](https://tomwenseleers.github.io/logconNMF/articles/algorithm.html#generalized-matrix-factorization)
relates its mean model to the GMF formulation.

**Deviance matrix factorization.** This is a rank-penalized, nonnegative
deviance matrix factorization on the mean scale, with hard log-concavity
constraints and an additive background. For standard exponential
families with fixed noise parameters, the loss is a weighted Bregman
divergence: Gaussian gives squared error, Poisson gives generalized KL,
and uncensored Gamma gives Itakura–Saito. See Wang and Carvalho (2023),
and Li and Tyagi (2023). The [mathematical
article](https://tomwenseleers.github.io/logconNMF/articles/algorithm.html#deviance-matrix-factorization-and-bregman-nmf)
gives the derivation and the negative-binomial case.

**Model and selection criterion.** With scans in rows, write
\$\$\mu=HS+B A+O,\qquad H,S,A\geq0.\$\$ Each profile has unit sum,
contiguous positive support and \\\Delta^2\log H\_{\cdot r}\leq0\\ on
that support. The objective is
\$\$J=-2\sum\_{ij}w\_{ij}\ell(y\_{ij};\mu\_{ij},\nu)+\lambda_R K,\$\$
where \\K\\ is the number of compounds and \\\nu\\ contains the shared
noise parameters. An IC preset uses \\\lambda_R=m\kappa\lambda\_{IC}\\,
with \\m\\ original channels. A numeric `lambda` supplies the total
compound penalty before multiplication by \\\kappa\\. The permissive
initialization uses the unadjusted preset \\\lambda\_{IC}\\; the
requested information factor is applied for final rank selection. Each
compound contributes one selection count; spectral entries and profile
complexity have no separate penalty. This is a compound-count
model-selection convention. The IC methods use \\mK\\ plus the number of
estimated shared noise parameters as their complexity count; it is
distinct from counting every free profile ordinate.

**Initialization and local updates.** A variance-matched grouped Poisson
working fit supplies candidate components, using
[`l0glmFit()`](https://tomwenseleers.github.io/L0glm/reference/l0glmFit.html)
with criterion-checked support refitting and compulsory background
columns. The subsequent updates use the specified family's exact score
and likelihood. For a spectral coefficient or profile ordinate, the
conditional mean is \\\mu_i=r_i+x_i b\\, with the other components fixed
in \\r_i\\. If \\g_i\\ and \\I_i\>0\\ are the negative log-likelihood
score and a positive working curvature, a nonnegative quadratic proposal
is \$\$b^+=\max\left\\0,b-\frac{\sum_i w_i x_i g_i} {\sum_i w_i x_i^2
I_i}\right\\.\$\$ Backtracking checks the actual loss, with an explicit
zero-boundary check. Poisson and quasipoisson use a specialized
bracketed scalar likelihood solve. Background terms use small
constrained quadratic solves, and all channels run independently through
OpenMP. Initial local profile proposals can use an adaptive ridge tether
centred on the current profile. Fine refinement instead uses constrained
Fisher scoring in log-profile coordinates, starting from the current
profile and spectrum. It enforces hard log concavity within the
quadratic update and checks the original-scale loss by backtracking.
This refinement adds no penalty toward zero or toward the initial PSF; a
small diagonal stabilizes the quadratic solve.

Spectral clustering proposes fusions without fitting a model. It
compares normalized spectra through Jensen–Shannon information loss and
joins overlapping compatible profiles. A least concave majorant of the
log profile enforces hard log concavity. Local refinement reallocates
signal between overlapping compounds and checks the actual criterion.
Every complete stage is accepted against its full working data, and
final refinement and scoring use all original observations.
`initialization` records the permissive starting search and its
unadjusted criterion. When the working-correlation factor changes the
final penalty, `history` starts from the resulting initialization and
records the final-criterion updates. Otherwise it also includes the
initialization stages. `ledger_error` checks the accumulated criterion
change within this ledger.

**Original-channel completion.** The final cleanup tests component
deletions, local fusions proposed by profile or spectral similarity, and
additional PSF candidates screened from original-channel residuals in
bounded scan windows. Deletion scores are cached and refreshed only for
overlapping profiles. Candidate construction performs no refitting.
Native local likelihood fits keep other compounds and the background in
the offset. Each accepted move decreases the complete original-channel
objective with unchanged numerical rank penalty, noise parameters,
observation weights and offsets. The `cleanup` entry records moves and
the complete objective change.

**Broad-profile fusion.** After refitting the original channels, an
additional algebraic pass compares elution profiles. Candidate pairs
must have cosine similarity at least 0.95 and at least one profile with
FWHM at least four times the typical PSF FWHM. This can reunite
complementary ion sets assigned to nearly proportional broad profiles.
The area-weighted profile is projected to hard log concavity and the
spectra are added. Each merge must improve the original-channel
likelihood-and-compound criterion with the noise scale, weights, offsets
and background fixed. Candidate fusion requires no refits. One final
local likelihood pass updates the profiles and spectra in affected
elution regions. It uses Fisher scoring in log-profile coordinates with
piecewise-linear interpolation and linear concavity constraints; the
noise likelihood is evaluated on the original intensity scale.
Backtracking preserves likelihood improvement. Neighboring profile modes
stay within one quarter of their starting FWHM, while merged profiles
can move freely. In log-profile coordinates \\\eta=\log H\_{\cdot r}\\,
the score and working information are \\g_t=H\_{tr}\sum_c
w\_{tc}S\_{rc}g\_{\mu,tc}\\ and \\I_t=H\_{tr}^2\sum_c
w\_{tc}S\_{rc}^2I\_{\mu,tc}\\. With interpolation matrix \\C\\, the
quadratic step minimizes \\\sum_t I_t\\(Cb)\_t-\eta_t+g_t/I_t\\^2\\,
subject to nonincreasing slopes of \\Cb\\. Zero-information rows receive
numerical stabilization. `profile_fusion` records the merge audit and
source component groups. Source IDs refer to components entering
original-channel profile refinement, or the final pass when that
refinement is disabled. Groups are retained across passes;
`source_components` follows the returned temporal ordering.

A second type of proposal absorbs a small, spectrally similar fragment
whose profile is contained within a broad component. Candidate spectra
have cosine similarity at least 0.97, at least 99\\ mass lies inside the
parent's 1\\ is at most 15\\ the representative PSF width. The broad
parent and a log-concave combined profile are tried as starts. Native
updates jointly refine each proposed component and its overlapping
neighbours, and the same refinement budget is applied to the unmerged
control. Acceptance uses the complete criterion on the original
channels, with noise settings, weights, offsets and the background
fixed. `contained_fusion` stores the trial diagnostics and source
groups.

**Alternative local starts.** Components with nearly coincident modes
can exchange signal during refinement. By default, final refinement also
tries doubling the left, right or both widths about the current mode,
interpolating the log profile to retain log concavity. Mode separation
must be at most one quarter of both starting FWHMs. Preliminary
screening and final acceptance both use every original channel to
preserve informative spectral contrasts. Promising starts are refined
locally; an unperturbed control receives at least the same refinement
budget. A restart is accepted only if it improves the complete
original-channel criterion relative to both the current fit and the
refined control, with component count, noise parameters and background
fixed. `profile_starts` records the proposals, accepted changes and
time. Component indices in its audit refer to the factors entering this
pass.

**Noise models.** Gaussian variance is \\\phi\\; quasipoisson variance
is \\\phi\mu\\; NB2 variance is \\\mu+\mu^2/\theta\\; standard Gamma
variance is \\\phi\mu^2\\. NB1 has variance \\\phi\mu\\, while
common-scale Gamma has shape \\\mu/b\\ and scale \\b\\, hence variance
\\b\mu\\. Their exact scores and densities are given in
[signal_families](https://tomwenseleers.github.io/L0glm/reference/signal_families.html).
A censored Gamma zero contributes a log CDF at the detection limit, and
a positive observation contributes its log density. For quasipoisson,
\\\ell_Q=\\y\log\mu-\mu-\log\Gamma(y+1)\\/\phi\\;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) marks this scaled
kernel with `quasi=TRUE`, and the regular IC method names return
quasi-information criteria. Compare candidate ranks under the same
family and frozen scale.

For an updated NB2 shape, the pooled score is
\$\$U\_\theta=\sum\_{ij}w\_{ij}\left\\\psi(y\_{ij}+\theta)-\psi(\theta)
+\log\theta+1-\log(\theta+\mu\_{ij})
-\frac{y\_{ij}+\theta}{\theta+\mu\_{ij}}\right\\.\$\$ Bounded Newton
proposals in \\\log\theta\\ are backtracked against the pooled
likelihood. One common shape is retained across original channels.

**Temporal calibration and aggregation.**
[`temporalDispersion()`](https://tomwenseleers.github.io/L0glm/reference/temporalDispersion.html)
uses fourth and eighth differences in quieter parts of a smoothed pilot
signal. This estimates the marginal noise scale; an optional lag-one
working correction multiplies the rank penalty by
\\\kappa=1+2\widehat\rho_1\\. It leaves the fitted variance parameter on
its marginal scale. This local covariance approximation is recorded in
`spec`; use `correlation="none"` for independent simulated data.

Preclustering uses
[`glmCluster()`](https://tomwenseleers.github.io/L0glm/reference/glmCluster.html)
with `size_fac=colSums(pmax(y,0))`, removing differences in total ion
intensity when comparing channel shapes. Poisson and common-parameter
NB1 latent signals add exactly; independent common-scale Gamma latent
signals also add. Standard NB2 and Gamma use moment-matched aggregate
parameters. Censored Gamma sums and differing channel weights give
working approximations in the reduced space. The returned model always
contains all original channels, their offsets and their original-channel
likelihood. The compound penalty uses the original channel count
throughout aggregation. Selection and initial fusion run in the reduced
space. Lifting the spectra to the original channels is followed by a
spectral/background refit and local shape refinement, with no repetition
of the reduced-space search. Worked full-channel fits to
[gcmsData1](https://tomwenseleers.github.io/logconNMF/reference/gcmsData1.md)
and
[gcmsData2](https://tomwenseleers.github.io/logconNMF/reference/gcmsData2.md),
with traces, heatmaps, profiles and spectra, are shown in the [GC/MS
deconvolution
tutorial](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
and the
[gcmsDeco()](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.html)
help page.
[writeGCMSCDF()](https://tomwenseleers.github.io/GCMSprime/reference/writeGCMSCDF.html)
exports these matrices for the complete file-based workflow.

**Smooth local profiles.** Original-channel refinement fits
\\\eta(t)=a+bt-\tfrac12\sum_j c_j\\(t-u_j)\_+^2-(t-v_j)\_+^2\\\\, with
\\c_j\geq0\\. Thus \\\eta''(t)\leq0\\ and \\H(t)=\exp\\\eta(t)\\\\ is
log-concave. Fisher scoring uses the original observation model and the
chain-rule score/information given above. Small constrained quadratic
solves and likelihood backtracking update each local profile. Joint
spectral updates separate overlapping compounds; compatible
broad-profile and tail fusions are accepted under the same
original-channel likelihood and compound-count penalty.

## References

Wang, L. and Carvalho, L. E. (2023). Deviance matrix factorization.
Electronic Journal of Statistics 17(2), 3762–3810.
[doi:10.1214/23-EJS2174](https://doi.org/10.1214/23-EJS2174) . Li, X.
and Tyagi, A. (2023). Block-active ADMM to minimize NMF with Bregman
divergences. Sensors 23(16), 7229.
[doi:10.3390/s23167229](https://doi.org/10.3390/s23167229) . Kidzinski,
L., Hui, F. K. C., Warton, D. I. and Hastie, T. J. (2022). Generalized
matrix factorization: efficient algorithms for fitting generalized
linear latent variable models to large data arrays. JMLR 23(291), 1–29.
<https://www.jmlr.org/papers/v23/20-1104.html>. Castiglione, C., Segers,
A., Clement, L. and Risso, D. (2026). Stochastic gradient descent
estimation of generalized matrix factorization models with application
to single-cell RNA sequencing data. Biostatistics 27(1), kxag010.
[doi:10.1093/biostatistics/kxag010](https://doi.org/10.1093/biostatistics/kxag010)
.

## See also

[gcmsDeco()](https://tomwenseleers.github.io/GCMSprime/reference/gcmsDeco.html),
[`l0glmBlindDeco()`](https://tomwenseleers.github.io/logconNMF/reference/l0glmBlindDeco.md),
[`l0glmFit()`](https://tomwenseleers.github.io/L0glm/reference/l0glmFit.html),
[`glmCluster()`](https://tomwenseleers.github.io/L0glm/reference/glmCluster.html)

## Examples

``` r
library(L0glm)
set.seed(12)
n <- 300L
psf <- convKernel("gaussian", 21L, list(sigma = 2))
X <- convMatrix(psf, n)
H <- cbind(X[, c(45, 95, 160, 235)],
           pmax(0, 1 - abs(seq_len(n) - 185) / 22))
H <- H %*% Matrix::Diagonal(x = 1 / Matrix::colSums(H))
S <- matrix(rexp(5 * 8, rate = 1 / 3000), 5, 8)
B <- splines::bs(seq_len(n), df = 4, intercept = TRUE)
mu <- as.matrix(H %*% S) + B %*% matrix(20, 4, 8)
Y <- matrix(rpois(length(mu), mu), n, 8)
fit <- splogconNMF(X, Y, background = B, family = poisson("identity"),
  precluster = 4L, correlation = "none", ncores = 2L)
fit
c(compounds = fit$rank, GIC = GIC(fit), seconds = fit$seconds[["total"]])
# Broad profiles with nearly identical shapes can share complementary ion sets.
# A conservative fusion pass runs on all original channels after preclustering.
fit$contained_fusion
matplot(fit$H, type = "l", lty = 1, xlab = "Scan", ylab = "Elution profile")

# For overdispersed intensities, estimate one common noise scale from temporal
# differences. Selection and fusion use four groups; the result has all channels.
Yq <- 5 * matrix(rpois(length(mu), mu / 5), n, 8)
fitq <- splogconNMF(X, Yq, background = B, family = quasipoisson("identity"),
  precluster = 4L, correlation = "none", ncores = 2L)
GIC(fitq)

# Other identity-link choices include gaussian(), mgcv::nb(), Gamma(),
# nb1(dispersion=5, link="identity"), and
# gammaScale(scale=5, link="identity", censor=0.5).

if (FALSE) { # \dontrun{
# Six co-eluting compounds, retaining all 671 mass channels.
data(gcmsData1)
Y <- gcmsData1
X <- convMatrix(convKernel("gaussian", 65L, list(sigma = 8)), nrow(Y))
f1 <- splogconNMF(X, Y, background = matrix(1, nrow(Y), 1L),
  family = quasipoisson("identity"), precluster = 16L,
  correlation = "auto", ncores = 8L)
c(compounds = f1$rank, GIC = GIC(f1), seconds = f1$seconds[["total"]])
matplot(as.numeric(rownames(Y)), as.matrix(f1$H), type = "l", lty = 1,
  xlab = "Retention time (min)", ylab = "Unit-area elution profile")

# Beer volatiles with a broad acid peak and a low-df spline background.
data(gcmsData2)
Y <- gcmsData2
time <- as.numeric(rownames(Y))
X <- convMatrix(convKernel("gaussian", 65L, list(sigma = 6)), nrow(Y))
B <- splines::bs(time, df = 4L, intercept = TRUE)
f2 <- splogconNMF(X, Y, background = B,
  family = quasipoisson("identity"), precluster = 32L,
  correlation = "auto", ncores = 8L)
c(compounds = f2$rank, GIC = GIC(f2), seconds = f2$seconds[["total"]])
plotMat(fitted(f2), xlab = "Retention time (min)", ylab = "Mass (m/z)",
  main = "Fitted GC/MS chromatogram")
# gcmsDeco() adds CDF input, MSP/Excel export, and traces/heatmaps/profile plots.
} # }
```
