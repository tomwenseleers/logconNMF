# Algorithms and theory for PSFs and sparse log-concave factorization

An ordered multichannel signal can be a mixture of overlapping
components. Each component has one profile over the ordered coordinate
and one vector of channel loadings. In GC–MS those factors describe how
the compound elutes and which ions it produces.

The algorithms below start with sparse deconvolution using a common
local response shape. Blind fitting or calibration can estimate that
response if needed. The log-concave factorization then allows each
retained component to have its own shape, and the simultaneous model can
share loadings across samples. Each stage uses the observation
likelihood and a criterion for model size.

This article follows the three tutorials: [PSF
estimation](https://tomwenseleers.github.io/logconNMF/articles/estimate-psf.md),
[single-sample
factorization](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.md)
and [simultaneous multisample
factorization](https://tomwenseleers.github.io/logconNMF/articles/simultaneous-nmf.md).
It specifies the actual objectives, coordinate conventions,
initialization and acceptance rules.

## Sparse spectra and log-concave profiles for count data

For a count matrix with time in rows and ions in columns, write the
expected signal in sample s as \mu\_{sij} = \epsilon + (B_s A_s)\_{ij} +
\sum\_{k=1}^K h\_{sik}p\_{kj}, \qquad p\_{kj}\ge0,\quad \sum_j p\_{kj}=1
. An absent compound has h\_{s\cdot k}=0. Otherwise its positive support
is a contiguous interval. For successive positive observations the
slopes (\log h\_{s,i+1,k}-\log h\_{sik})/(t\_{s,i+1}-t\_{si}) are
nonincreasing. This permits asymmetric peaks, including an abrupt change
of slope at the apex. Each h\_{sik} is the compound’s expected total ion
count in that scan. Two samples can have different peak shapes, widths,
areas and acquisition times.

With \operatorname{Var}(Y\_{sij})=\phi_s\mu\_{sij}, the conditional
updates minimize J =
2\sum\_{sij}\frac{y\_{sij}\log(y\_{sij}/\mu\_{sij})-y\_{sij}+\mu\_{sij}}{\phi_s} +
\lambda\_{\rm eff} \left\[\sum\_{k=1}^{K}(\lVert
p_k\rVert_0-1)+\sum\_{sk}d(h\_{s\cdot k})\right\]. Use 0\log 0=0. A
spectrum with u positive ions contributes u-1 free proportions. A
positive log profile on two or more scans contributes two endpoint
ordinates plus one dimension per interior strict decrease in log slope.
A single-scan support contributes one dimension and an absent profile
contributes zero. This is the dimension of the fitted face of the
log-concavity cone, conditional on support and knots. Searching over
their locations adds selection uncertainty. The working GIC uses
\lambda\_{\rm eff}=\kappa\log(p)\log\log(n). For a joint fit use
p=n=\sum_s n_s. The temporal information factor \kappa equals one for
independent residuals. Backgrounds are estimated on a fixed basis and
remain unpenalized. The fitted noise parameters add a constant charge to
the reported information criterion for that dataset.

1.  **Form a spectral dictionary.** Calibrate the PSF, form its banded
    design and fit a nonnegative grouped L0 regression across channels.
    For a single sample the independent log-concave fit supplies local
    starts. Repeated acquisitions can share a PSF estimated from pooled,
    count-conserving total-ion traces. Channel aggregation accelerates
    the initial regression. Estimate the original-ion intensities by
    nonnegative least squares. The sparse Gram matrix records which
    initial profiles overlap.
2.  **Locate candidate occurrences.** With several samples, interpolate
    each peak’s 1%-height onset on its native scan grid. Reciprocal
    spectral matches within the retention window provide the anchors for
    monotone time maps. On that shared clock, propose local spectral
    directions that overlapping candidates cannot already explain by a
    nonnegative mixture. Retain native time coordinates for testing peak
    shape. With one sample the time map is the identity.
3.  **Alternate two conditional updates.** Visit components in
    decreasing total abundance. To update one profile, put all other
    components and background in the offset. Compute per-scan
    conditional count estimates, project their logarithms onto the
    concavity cone with information weights, optimize the amplitude, and
    compare with keeping or removing the old profile using original-ion
    J. A curve reaching a window boundary can expand. For fixed
    profiles, refit a common spectrum from all occurrences. Compare each
    ion with zero using the same penalty and transfer the spectrum’s
    normalization into the profiles. Update the background.
4.  **Check ambiguous overlaps.** Test proposed local support exchanges
    and joint refits against the same original-channel objective. The
    default allows 12 initial sweeps and up to two exchange rounds, with
    at most ten refit sweeps after an accepted exchange. Relative
    objective improvement below 10^{-3} ends a round. These bounded
    searches return a local solution. Supplied ground truth need not
    minimize a penalized criterion on one realization of noisy data.

For overlapping peaks, a spectrum observed in relatively few samples can
be proposed as a nonnegative combination of recurring spectra from the
same aligned retention interval. Each proposed template must also occur
in a sample where the putative mixture is absent. The solver refits the
contributing log-concave profiles and shared spectra on the original ion
counts. It tests both redistribution and retention of a residual
spectrum. Within a sample, profiles whose half-height intervals overlap
and whose spectra are nearly linearly dependent are also refitted
jointly. The fit compares keeping the locally redundant profiles with
dropping them. The fitted background and other compounds contribute to
the offset in each comparison. An exchange is accepted only when it
lowers the same deviance-plus-dimension penalty as the factor updates.
Library identities and known peak areas enter only at validation.

With the spectra fixed, profiles in different samples have separate
likelihood contributions. With profiles fixed, the spectral likelihood
and ion inclusion penalty separate by ion. Workers therefore update
different samples or ion columns while preserving the component order
within each. The spectral sweep transfers the unit-sum normalization
into each profile after all ions have finished. Accepting a merge
changes the offsets of its neighbours, so merge decisions remain
sequential.

For simulated data a truth-based comparison has a precise
interpretation. Write \mu^\star=E(Y) and fix dispersion and weights. The
expected difference in Poisson deviances between any fixed \mu and
\mu^\star is 2\sum\[\mu^\star\log(\mu^\star/\mu)+\mu-\mu^\star\]/\phi,
which is nonnegative. On a single noisy dataset this inequality need not
hold. Adding a dimension penalty can also favor a simpler approximation
in expectation. A true-profile initialization consequently provides a
reference value for this dataset. A certificate of global optimality
would require a lower bound over all admissible factorizations. In a
grouped Poisson regression, starting the fitted mean at HS+BA gives the
correct first Fisher weights. The support penalty and the sample score
can still move the coefficients away from their generating values.

Use `objective.trace=TRUE` to record original-channel deviance, the two
dimension counts, rank, occurrences and J after each block. The values
in `conditional_objective` use the saturated model as the deviance
reference. The
[`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
extractor includes the data-dependent likelihood constant and estimated
noise dimensions. Comparisons on the same observations and noise
specification give the same objective differences.

### Convergence from estimated and true elution profiles

On the small 16-sample benchmark (1,200 scans and 258 ions per sample),
the standard fit recovered all 769 true occurrences with only 26 false
positives. After dictionary construction its GIC fell from 21,188,577 to
4,173,728. Of the total decrease, 91.9% occurred in the first three
profile and spectrum sweeps and 96.7% in the first five. Across 56
checkpoints there were 51 decreases and four ties.

Starting the same joint updates from the true unit-area profiles and
their known cross-sample matches gave a final GIC of 4,154,879. The
standard fit retained 115 shared factors and the reference fit retained
100. The standard fit finished **0.454% above this empirical
reference**. Scaled deviance accounted for 80% of the 18,848 GIC units
between the endpoints and model complexity for the remaining 20%. The
signed percentage is 100(\mathrm{GIC}\_{\rm standard}-\mathrm{GIC}\_{\rm
reference})/ \|\mathrm{GIC}\_{\rm reference}\|. Either start can give
the lower value. Certifying a global optimum for this nonconvex problem
would require a lower bound on the criterion.

![GIC over the joint fit and percentage difference from the
truth-initialized endpoint on the small 16-sample
benchmark](../reference/figures/empirical-gic-convergence.png)

GIC over the joint fit and percentage difference from the
truth-initialized endpoint on the small 16-sample benchmark

[Download the block-by-block GIC, deviance, rank and degrees of freedom
for both
fits](https://tomwenseleers.github.io/logconNMF/benchmark-gic-trace.csv).

For a simulation from
[`simulateGCMSSamples()`](https://tomwenseleers.github.io/logconNMF/reference/simulateGCMSSamples.md),
pass the simulation itself or a named list of scan-by-compound matrices
with shared compound IDs. The diagnostic estimates initial areas,
spectra and backgrounds from the counts and then applies the same
profile, spectral and support updates to both starts. Irregular sampling
can make the true integrated counts slightly nonconcave on the observed
time grid. The reference takes the least concave majorant of log height
on each support interval and rescales it to unit area. On this benchmark
the median fraction of a profile’s area redistributed by this operation
was 0.002%, and the largest was 0.044%.

``` r

library(logconNMF)
sim <- simulateGCMSSamples(samples = 2L, n = 240L, k = 8L, seed = 3L)
samples <- lapply(sim$samples, function(s)
  list(y = s$Y, time = s$time, sample_id = s$sample_id))
fit <- logconNMF(samples, family = quasipoisson("identity"), truth.profiles = sim)
gicGap(fit)
plot(fit$gic_gap)
```

With profile shapes held fixed, a single spectral sweep on this
simulation already achieved a median cosine of **0.99999895** against
the 100 true spectra. That sweep also updates the background and
transfers each spectrum’s total scale to its profiles. The full
reference fit additionally updates each sample’s shape and area. Its
scaled deviance reached 3.523 million, compared with 5.833 million after
the first fixed-shape spectral sweep and 5.732 million after 16
fixed-shape sweeps. Overlapping compounds share observations, so
updating one spectrum changes the residual available to its neighbours.

If the complete profile matrices H_s (including areas) and backgrounds
are fixed, estimating the shared spectrum matrix P is one conditional
regression problem. With fixed support, the Poisson Hessian for ion j is
\sum_s
H_s^\mathsf{T}\operatorname{diag}(y\_{sij}/(\phi_s\mu\_{sij}^2))H_s,
which is positive semidefinite. A numerical solver can take several
sweeps to solve that single block. Supplying only unit-area shapes also
leaves a separate amplitude to estimate for each sample and compound.

With eight threads, the complete CDF workflow averaged 50.9 s with the
diagnostic off, 57.0 s with only tracing, and 76.2 s with the reference
fit and both traces (two runs per setting). The extra fit reuses the
estimated PSF, registration and noise calibration. Ordinary fits leave
tracing and the reference calculation off.

The conditional solver supports identity-link Poisson and quasipoisson
data with unit observation weights and exposure. Its C++ kernels update
sample profiles and independent ion columns with OpenMP. Fits use the
arrays held by R. For on-disk inputs, the solver reads the stored counts
for each sample. The general solver supplies the other identity-link
families, observation weights, offsets and exposures.

## Rank selection in the general identity-link solver

For a fitted mean M_s in each sample, `algorithm="classic"` searches for
nonnegative channel loadings and nonnegative, log-concave profiles that
approximately minimize

\begin{aligned} \mathcal J &=
-2\sum_s\ell\_{\nu_s}(M_s;Y_s)+\kappa\lambda d\_{\mathrm{sel}},\\
d\_{\mathrm{sel}} &= d\_\nu + \begin{cases} mK, & \text{one sample},\\
(m-1)G+L, & \text{simultaneous samples}. \end{cases} \end{aligned}

Here m is the original channel count, K is the single-sample rank, G
counts shared spectra and L counts detected sample–component
occurrences. A normalized shared spectrum has m-1 free coordinates. Its
amplitude in each sample contributes one further coordinate, so G=L=K
gives (m-1)K+K=mK. This is a selection penalty on components and
occurrences. Profile shapes and smooth backgrounds have no separate
charge. The fit can use clustered channels for initialization. The
likelihood in the final comparison uses all original channels.

The weighted likelihood and quasi-Poisson score are defined below. The
optional working temporal adjustment multiplies the penalty by \kappa,
which equals one under independence. The default penalty
\lambda=\log(p)\log\log(n) uses the dictionary width and scan count for
one sample. For simultaneous fits, both preset inputs equal the total
scan count, p=n=\sum_s n_s. A numeric `lambda` uses the same units in
both modes. The reported criteria and simultaneous objective count
estimated noise parameters through d\_\nu. The single-sample `objective`
omits this additive constant for a fixed noise specification.

## Observation model and notation

Let Y=(y\_{ij})\in\mathbb R^{n\times m} contain ordered observations,
with scans or time points in rows and response channels in columns. Let
M=(\mu\_{ij}) be the nonnegative mean signal and \nu the shared noise
parameters. We use the weighted log likelihood

\ell\_\nu(M;Y)=\sum\_{i=1}^{n}\sum\_{j=1}^{m}
w\_{ij}\ell\_\nu(y\_{ij};\mu\_{ij}),\qquad w\_{ij}\geq0.

For quasipoisson, \ell\_\nu denotes the scaled Poisson quasi-likelihood
kernel defined in the [noise-model section](#noise-models). Observation
weights are applied once in this definition. The seven supported
identity-link models use their own mean-variance relationships and
likelihood derivatives; estimation and acceptance take place on the
original observation scale.

Throughout, n is the number of original observations, m the original
channel count and p the number of candidate PSF columns. B is a
nonnegative low-dimensional background basis, A\geq0 its coefficients,
and O known additive signal. The notation matches the [GCMSprime
mathematical
basis](https://tomwenseleers.github.io/GCMSprime/articles/algorithm.html),
which develops the chromatographic interpretation of these models.

## Convolution and PSF estimation

### Deconvolution based on a known peak shape

Let C\in\mathbb R\_+^{p\times m} be a spike-amplitude matrix and g a
known point spread function (PSF), or peak shape function. The
fixed-kernel mean model is

M=X(g)C+BA+O,\qquad C,A\geq0.

For a constant kernel, X(g)\_{ir}=g\_{i-r}, with the specified centering
and boundary convention. A nonzero row calls a candidate peak using
information from every response channel. The grouped L_0 selection
target is

\begin{aligned} s(C)&=\\C\\\_{2,0}
=\sum\_{r=1}^{p}\mathbf1\\\\C\_{r\cdot}\\\_2\>0\\,\\ \mathcal G(C,A;g)
&=-2\ell\_\nu\\X(g)C+BA+O;Y\\+\lambda_g s(C). \end{aligned}

Here \lambda_g is the total row-selection charge. This is the target for
**best subset of groups selection (BSGS)**, where a group consists of
all channel coefficients at one candidate position. Independent
selection instead uses \sum\_{rj}\mathbf1\\C\_{rj}\ne0\\. Background
terms are compulsory and unpenalized in the coefficient updates.

L0glm approximates these sparse-selection objectives with adaptive ridge
within GLM IRLS, followed by optional relaxed refitting. Its GIC preset
supplies the information-criterion scale

\lambda\_{\rm IC}=\log p\\\log\log n,

with natural logarithms; AIC uses 2 and BIC uses \log n. The
high-dimensional GIC scale follows [Fan & Tang
(2013)](https://doi.org/10.1111/rssb.12001). The working adaptive ridge
uses a quarter of the selection scale with the appropriate dispersion
factor. See [L0glm’s
algorithm](https://tomwenseleers.github.io/L0glm/articles/algorithm.html)
for the inverse-square updates and their quadratic calibration. The
grouped candidate fit supplies a support for the factorization; its
final [rank-selection criterion](#rank-selection) has an explicit
original-channel compound charge.

The sparse pulse-train formulation was inspired by [de Rooi & Eilers
(2011)](https://doi.org/10.1016/j.aca.2011.05.030). Our extensions add
the information-criterion regularization preset, grouped selection
approximating BSGS, and Poisson or overdispersed Poisson observation
models. With estimated quasipoisson dispersion, the corresponding target
is quasi-GIC, using the same row-support penalty and a scaled
quasi-likelihood.

### Relation to convolutive NMF

[Cichocki et al.](https://doi.org/10.1002/9780470747278) write
convolutive NMF (CNMF) as

Y=\sum\_{\ell=0}^{L-1}A\_\ell Z^{(\ell)}+E,

where Z^{(\ell)} shifts a nonnegative source matrix along its ordered
coordinate and A\_\ell\geq0 is a lag-specific basis matrix. The
constant-PSF model above has the exact transposed representation

\\X(g)C\\^{\mathsf T} =\sum\_\ell (g\_\ell I_m)(C^{\mathsf T})^{(\ell)}.

It is therefore a constrained CNMF with A\_\ell=g\_\ell I_m: all
channels share one scalar kernel and there is no lag-dependent channel
mixing. If the spike-amplitude matrix is further factorized as C=ZP,
where P contains component spectra, then

\\X(g)ZP\\^{\mathsf T} =\sum\_\ell \underbrace{(g\_\ell P^{\mathsf
T})}\_{A\_\ell} (Z^{\mathsf T})^{(\ell)}.

Thus blind constant-PSF estimation is a highly constrained,
lag-separable convolutive-NMF problem as well as a sparse multichannel
blind-deconvolution problem. This classification does not extend
indiscriminately to every stage: calibration from isolated standards
estimates a kernel without an NMF; smoothly time-varying PSFs are
shift-variant and only locally convolutional; and the final
factorization estimates independent component- and sample-specific
log-concave profiles rather than shifted copies of one kernel. Here
**CNMF** is reserved for *convolutive NMF*.

### Calibration with isolated peaks

For a local isolated standard peak, write

\mu\_{ij}=b\_{ij}+o\_{ij}+h_i s_j, \qquad h_i\>0,\quad s_j\geq0,\quad
D^2\log h\leq0.

A low-quantile spline estimates background. Signal-to-noise, separation
and profile geometry identify informative peak windows. A native
constrained log-profile update fits each window on the original
observation scale. Normalize its profile, align by centroid, and
resample onto the sample’s median scan interval. A high quantile of
fitted profile support widths supplies the common support when
`control$L` is omitted. `control$L` fixes an odd length.

**Algorithm P1 — standard PSF calibration**

1.  Resolve the shared noise parameters, weights and additive offsets.
    Optionally aggregate channels; carry the corresponding sum-model
    parameters.
2.  Estimate a lower-quantile background and detect informative isolated
    peaks.
3.  For each retained local window, initialize s_j from positive
    integrated signal and h from its positive total trace projected to
    log concavity.
4.  Update \eta=\log h under D^2\eta\leq0, using the likelihood score
    and positive working curvature. Backtrack against the original
    likelihood.
5.  Set h\leftarrow h/\sum_i h_i, align its centroid, and place it on
    the common lag grid. Infer support from fitted peak mass or use the
    supplied length.
6.  Interpolate the calibrated log kernels at the requested PSF anchor
    positions, normalize and construct the nonnegative sparse
    convolution design.

### Alternating blind estimation

For a fixed channel spike train c_j, convolution is linear in the
kernel:

X(g)c_j=T(c_j)g,

where T(c_j) contains shifted copies of the spike train. This gives a
second nonnegative GLM. Blind estimation therefore alternates **two
L0glm fits**: one estimates spike trains and background given the
current kernel; the other estimates the kernel given those spike trains.
They form an alternating block-coordinate procedure for the joint
sparse-likelihood target

\mathcal G(C,A,g) =-2\ell\_\nu\\X(g)C+BA+O;Y\\+\lambda_g s(C),

under nonnegativity and a normalized, log-concave kernel. The
conditional updates use adaptive ridge and kernel shape smoothing. Their
complete spike/background reconstructions are compared using the
requested IC.

The PSF fit’s ridge center is the current kernel. The implementation
normalizes each informative channel’s spike train for conditioning and
rescales the prior center consistently. Positive starts retain kernel
tails. Its default kernel update keeps alternating lag positions as
compulsory terms; the remaining positions use adaptive selection. A
log-concave smoother and centroid alignment produce the kernel proposal,
which is normalized to unit sum.

Completed proposals are compared after re-estimating their spike trains
and backgrounds. A full-tail proposal and a pedestal-removed proposal
can both be considered. The returned checkpoint is the synchronized
spike/background fit with the lowest requested criterion among these
completed candidates.

Both fits use the same observation model and resolved noise parameters.
For fixed-dispersion quasipoisson counts, weighting the Poisson log
likelihood by 1/\phi gives the same mean-dependent quasi-likelihood.
This allows the native Poisson backend to retain the calibrated noise
scale. The other families use their corresponding likelihoods. The
kernel updates and completed checkpoint comparisons approximate the
joint selection problem.

**Algorithm P2 — blind constant-PSF estimation**

1.  Set g^{(0)} to the supplied/calibrated kernel, or a normalized
    Gaussian on the requested lag grid. Resolve shared noise parameters
    once in
    [`estimatePSF()`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSF.md).
2.  For t=0,\ldots,T-1:
    1.  Construct X(g^{(t)}). Fit \[X(g^{(t)}),B\] to Y, with compulsory
        background and the chosen sparse/grouped GLM preset.
    2.  Store this spike/background fit and its criterion.
    3.  For every informative channel, construct T(c_j) and fit its
        kernel using background plus known signal as offset, the current
        kernel as ridge center, and the same observation weights/noise
        model.
    4.  Smooth, recenter and normalize the channel kernels; combine
        informative channel estimates. Construct the optional
        concentrated alternative.
    5.  Refit spike/background models for the competing proposals and
        select the better completed candidate.
    6.  Stop if \\g^{(t+1)}-g^{(t)}\\\_1/\\g^{(t)}\\\_1\<\epsilon.
3.  Synchronize the last kernel with a spike/background refit. Return
    the best scored checkpoint, its sparse design and its selection
    trace.

`control$maxiter` bounds this loop and `control$tol` is \epsilon.
`method="auto"` runs P1 followed by P2; `"standard"` runs P1 and
`"blind"` runs P2. A supplied kernel with `maxiter=0L` constructs its
design directly.

### Smoothly time-varying kernels

With `nwindows>1L`, run local estimation in overlapping windows centered
at \tau_1\<\cdots\<\tau_R. Use a common lag grid and interpolate
adjacent log kernels:

\eta_t=(1-a)\log g_r+a\log g\_{r+1},\quad a=u^2(3-2u),\quad
u=\frac{t-\tau_r}{\tau\_{r+1}-\tau_r}, \qquad
g_t=\frac{\exp\eta_t}{\sum\_\ell\exp\eta\_{t\ell}}.

A convex combination of concave log profiles is concave; normalization
changes only its intercept. The smoothstep weights give zero
interpolation slope at anchor joins. Each design column uses its
position’s kernel, is projected to contiguous log-concave support, and
is normalized after boundary truncation. The measured observations stay
on their scan grid; kernels are indexed in scans.

**Algorithm P3 — windowed PSF estimation**

1.  Choose evenly spaced sample anchor positions and overlapping local
    windows.
2.  Run P1/P2 for each anchor, retaining a common lag grid and its
    diagnostics.
3.  Interpolate log kernels using the weights above, with endpoint
    continuation.
4.  Place each local kernel at its candidate impulse position, truncate
    at the observation boundary, enforce log concavity and normalize the
    column.
5.  Return sparse X, anchor kernels, times, widths and convergence
    diagnostics.

These PSFs describe a common response at each position. The
factorization below subsequently estimates component-specific shapes.

## Single-sample model and rank criterion

For scans i, channels j and compounds r,

\mu\_{ij}=\sum\_{r=1}^K H\_{ir}S\_{rj}+(BA)\_{ij}+O\_{ij},\qquad
H,S,A\geq0.

The columns of H have unit sum, contiguous support and \Delta^2\log
H\_{\cdot r}\leq0. Spectra S therefore have the interpretation of
discrete integrated channel intensities. The objective is

J(H,S,A,K)=-2\ell\_\nu(HS+BA+O;Y)+\lambda_RK, \qquad
\lambda_R=m\kappa\lambda\_{\rm IC}.

Here m is the original channel count, \nu contains pooled noise
parameters and \kappa=1 for independent observations. The GIC preset
uses \lambda\_{\rm IC}=\log p\\\log\log n, resolved from the original
PSF design and observation count. For numeric `lambda` equal to c in
[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md),
the compound charge is \lambda_R=m\kappa c. The single-sample solver
[`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
accepts the total compound charge before temporal adjustment.

Dividing by the original channel count gives the equivalent criterion

J/m=-2\ell\_\nu/m+\kappa\lambda\_{\rm IC}K.

Thus each compound contributes one rank-selection unit on the
average-channel likelihood scale. Background coefficients and profile
complexity receive no additional penalty under this compound-count
convention. Channel aggregation retains the original m in \lambda_R;
changing the number of outcome groups does not redefine the compound
penalty.

### Reported information criteria

Let d\_\nu count estimated shared noise parameters. The reported IC
methods use the selection complexity d\_{\rm sel}=mK+d\_\nu:

\operatorname{IC}=-2\ell\_\nu +\kappa\lambda\_{\rm IC}(mK+d\_\nu).

For a fit using the matching IC preset, this equals
J+\kappa\lambda\_{\rm IC}d\_\nu. The extra term is constant when
comparing ranks with the same noise specification. Known noise
parameters add no selection count. The reported count implements the
compound-selection convention; it leaves profile ordinates and smooth
background coefficients unpenalized.
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html),
[`EBIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
and
[`MBIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
use the fitted likelihood and this count. Quasipoisson returns the
corresponding quasi-information criteria through the same method names.

### Log-concave component-specific profiles

On a regular grid, the shape constraint is H\_{ir}^2\geq
H\_{i-1,r}H\_{i+1,r} at interior positive ordinates. On an irregular
grid, successive slopes of \log H\_{\cdot r} in actual time must be
nonincreasing. Each component has contiguous support. This accommodates
Gaussian, exponentially modified Gaussian and asymmetric triangular or
shark-fin-shaped chromatographic profiles. Gaussian and exponential
profiles are log-concave, and convolution preserves this property
([Saumard & Wellner, 2014](https://doi.org/10.1214/14-SS107)).

Log-concavity implies unimodality and supplies additional control over
the rising and falling slopes. It provides stronger shape information
than the unimodality constraints commonly used in MCR-ALS ([de Juan,
Jaumot & Tauler, 2014](https://doi.org/10.1039/C4AY00571F)). Estimating
each profile separately allows compound- and concentration-dependent
broadening and tailing around the common PSF-based initialization.

## Generalized matrix factorization

The broad statistical setting is **generalized matrix factorization
(GMF)**: an exponential-family mean is linked to a low-rank product,
optionally with observed covariate effects and offsets. See [Kidziński
et al. (2022)](https://www.jmlr.org/papers/v23/20-1104.html) and
[Castiglione et
al. (2026)](https://doi.org/10.1093/biostatistics/kxag010). For example,
the latter mean-model structure is

g(\mu)=X\_{\rm row}B\_{\rm row}^{\top} +B\_{\rm col}X\_{\rm
col}^{\top}+UV^{\top}+O.

Our model uses the identity link, U=H, V^{\top}=S, and X\_{\rm row}=B,
B\_{\rm row}^{\top}=A for the smooth background; the column-covariate
term is omitted. It adds nonnegativity, unit-area log-concave profiles
and a compound-count penalty to select the rank. This gives a sparse
log-concave nonnegative GMF with automatic rank selection.

## Deviance matrix factorization and Bregman NMF

[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
fits a **rank-penalized, nonnegative deviance matrix factorization with
log-concave profiles**. It belongs to the deviance-based framework of
[Wang & Carvalho (2023)](https://doi.org/10.1214/23-EJS2174), using an
identity link for the additive mean \mu=HS+BA+O. The loss also connects
directly to [NMF with Bregman
divergences](https://doi.org/10.3390/s23167229).

For a regular exponential-dispersion model with fixed noise parameters,

\ell(y;\vartheta,\phi)=\frac{y\vartheta-b(\vartheta)}{\phi}+c(y,\phi),
\qquad \mu=b'(\vartheta),

let \Psi=b^\* be the convex conjugate of b. Its Bregman divergence is

D\_\Psi(y,\mu)=\Psi(y)-\Psi(\mu)-\Psi'(\mu)(y-\mu).

Since \Psi'(\mu)=\vartheta(\mu) and
\Psi(\mu)=\mu\vartheta(\mu)-b\\\vartheta(\mu)\\, substitution gives

\begin{aligned} D\_\Psi(y,\mu)
&=y\vartheta(y)-b\\\vartheta(y)\\-y\vartheta(\mu)+b\\\vartheta(\mu)\\,\\
2\\\ell(y;y,\phi)-\ell(y;\mu,\phi)\\ &=\frac{2}{\phi}D\_\Psi(y,\mu).
\end{aligned}

Here the likelihood on the second line is written in mean coordinates.
Boundary observations use the corresponding limits. Thus, conditional on
the noise parameters, fitting is equivalent to minimizing

2\sum\_{ij}\frac{w\_{ij}}{\phi}D\_\Psi(y\_{ij},(HS+BA+O)\_{ij})
+\lambda_RK

under the nonnegativity, profile-normalization and log-concavity
constraints. The saturated likelihood contributes a constant in these
conditional model comparisons. The corresponding divergences include:

| Noise model | Generator \Psi(u) | Half unit deviance D\_\Psi(y,\mu) |
|:---|:---|:---|
| Gaussian | u^2/2 | (y-\mu)^2/2 |
| Poisson | u\log u-u | y\log(y/\mu)-y+\mu |
| Gamma, uncensored y\>0 | -\log u | y/\mu-\log(y/\mu)-1 |
| NB2, fixed shape \theta | u\log u-(u+\theta)\log(u+\theta) | y\log(y/\mu)-(y+\theta)\log\\(y+\theta)/(\mu+\theta)\\ |

For NB2, \phi=1 and the shape \theta enters the generator. Gaussian,
Poisson and Gamma therefore give squared-error, generalized KL and
Itakura–Saito NMF losses. Quasipoisson gives the KL kernel scaled by its
dispersion. Censored Gamma and the mean-dependent-shape extensions use
their full likelihood kernels described below; estimated noise
parameters are compared using the full likelihood as well.

## Smooth background and flexible noise models

The additive background BA uses a nonnegative, low-dimensional B-spline
basis, or a constant basis for a short window. Its coefficients are
compulsory and unpenalized. Conditional background updates retain the
fitted component signal in the offset and use the specified observation
likelihood. For chromatograms the basis should vary much more slowly
than an individual peak. A lower-quantile background fit initializes
isolated-peak calibration.

### Seven identity-link signal models

| Noise model | Variance | Negative log-likelihood score in the mean |
|:---|:---|:---|
| Gaussian | \phi | (\mu-y)/\phi |
| Poisson | \mu | 1-y/\mu |
| Quasipoisson | \phi\mu | (1-y/\mu)/\phi |
| NB2 | \mu+\mu^2/\theta | \theta(\mu-y)/\\\mu(\mu+\theta)\\ |
| Gamma | \phi\mu^2 | (\mu-y)/(\phi\mu^2) |
| NB1 | \phi\mu | \\\psi(a)-\psi(a+y)+\log\phi\\/(\phi-1), a=\mu/(\phi-1) |
| Common-scale Gamma | b_0\mu | \\\psi(\mu/b_0)+\log b_0-\log y\\/b_0 |

The final two models use the exact distributions documented in
[`signal_families`](https://tomwenseleers.github.io/L0glm/reference/signal_families.html).
Their observed positive curvatures are differences of trigamma functions
for NB1 and \psi_1(\mu/b_0)/b_0^2 for common-scale Gamma. The latter is
convex in the mean; the NB1 curvature is nonnegative because trigamma
decreases on the positive real line. At \phi=1, NB1 uses its Poisson
limit. The Gamma scores in the table apply to positive, uncensored
observations; zero observations use the censoring likelihood below.
Gaussian variance is \sigma^2=\phi; common-scale Gamma uses b_0 for its
shared scale.

Poisson, NB1 and NB2 require integer counts. Quasipoisson accepts
nonnegative continuous intensities; Gaussian observations may also be
negative, while the modeled component signal remains nonnegative. Gamma
observations below the detection threshold are represented by zero.

[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
defaults to `poisson("identity")` for count-data factorization.
[`estimatePSF()`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSF.md)
defaults to `quasipoisson("identity")`. For GC/MS electron-multiplier
intensities, GCMSprime uses the latter as its workflow default and
reuses the estimated noise specification during PSF estimation and
factorization.

### Quasipoisson and estimated overdispersion

Poisson ion arrivals combined with random electron-multiplier gain
motivate an overdispersed counting model. For independent identically
distributed gains G_a, independent of N,

Y=\sum\_{a=1}^{N}G_a,\quad N\sim\operatorname{Poisson}(\Lambda),\qquad
\operatorname{E}(Y)=\Lambda\operatorname{E}(G),\quad
\operatorname{Var}(Y)=\Lambda\operatorname{E}(G^2).

The variance-to-mean ratio depends on gain and intensity units. A
quasipoisson working model retains this proportional variance structure;
electronics and digitization can add further variation. See [Ipsen
(2015)](https://doi.org/10.1021/ac503554u) for the physical measurement
model.

With a resolved common dispersion \widehat\phi, the implemented scaled
kernel is

\ell_Q(M;Y)=\frac{1}{\widehat\phi}\sum\_{ij}w\_{ij}
\\y\_{ij}\log\mu\_{ij}-\mu\_{ij}-\log\Gamma(y\_{ij}+1)\\.

The final term is constant for fixed observations and scale. It
therefore gives the same model comparisons as the mean-dependent
quasi-likelihood in GCMSprime. Equivalently, write

\begin{aligned} D_P(Y,M)&=2\sum\_{ij}w\_{ij}
\\y\_{ij}\log(y\_{ij}/\mu\_{ij})-y\_{ij}+\mu\_{ij}\\,\\
J_Q&=D_P(Y,M)/\widehat\phi+\lambda_RK+\text{constant}. \end{aligned}

Use 0\log(0/\mu)=0. Thus the same original-channel rank charge applies
to GIC and quasi-GIC; the estimated dispersion scales the fit term.
[`logLik()`](https://rdrr.io/r/stats/logLik.html) marks the quasipoisson
kernel with `quasi=TRUE`, and the regular IC method names return
quasi-information criteria.

### Pooled noise calibration

Noise parameters are shared across channels. With `dispersion=NULL`, the
marginal variance or scale is estimated from fourth- and eighth-order
temporal differences in quieter regions of a smoothed pilot signal.
Supplying a positive dispersion fixes it. NB2 similarly supports a fixed
positive `theta`, automatic pooled initialization with `theta=NULL`, or
a negative starting value that enables subsequent likelihood-improving
shape updates. Poisson and NB2 have GLM scale dispersion one; NB2
overdispersion is represented by \theta.

Candidate rank or grouping moves are compared using the same resolved
noise parameters, observation weights, offsets and numerical penalty. An
enabled NB2 shape update is a separate likelihood-improving step; the
next candidate comparison uses that updated common shape. This keeps
each acceptance test on a common noise scale. The [temporal-variation
section](#temporal-variation) details the marginal variance and optional
correlation adjustment.

### Censored Gamma intensities

Recorded zeros may represent intensities below a detection limit d,
whose contribution is \log P(Y\<d;\mu). Positive observations contribute
their ordinary density. For standard Gamma, set a=1/\phi, z=d/(\phi\mu);
its negative log-CDF derivative is

g(\mu)=\frac{z^a e^{-z}}{\Gamma(a)P(a,z)\mu}.

For common-scale Gamma, z=d/b_0 stays fixed and a=\mu/b_0 changes. The
convergent incomplete-Gamma series

P(a,z)=\frac{e^{-z}z^a}{\Gamma(a+1)}
\left\\1+\sum\_{k\geq1}\frac{z^k}{(a+1)\cdots(a+k)}\right\\

provides its shape derivative when z is small. The general branch uses
the Gamma CDF and a numerical shape derivative. Both branches are
checked against the actual censored likelihood. The default detection
limit is 0.5 intensity units and can be changed through
`control$censor`.

## Local quadratic updates

For a coefficient b\geq0, keep the rest of the signal in an offset r_i
and write \mu_i=r_i+x_i b. With negative log-likelihood score g_i and
positive working curvature I_i, propose

b^+=\max\left\\0,b-\frac{\sum_i w_i x_i g_i} {\sum_i w_i x_i^2
I_i}\right\\.

Actual-loss backtracking and a zero-boundary check complete each scalar
solve. The Poisson branch uses a faster bracketed likelihood solve. Its
score \sum_iw_ix_i(1-y_i/\mu_i) is monotone, and its nonnegative root is
bounded above by \sum_iw_iy_i/\sum_iw_ix_i. Channels are independent
conditional on the profiles, allowing deterministic OpenMP parallelism.

Local profile updates use the original observation likelihood. For
\eta_i=\log h_i, \mu\_{ij}=r\_{ij}+h_i s_j, the chain-rule score and
working information are

g\_{\eta,i}=h_i\sum_jw\_{ij}s_jg\_{\mu,ij},\qquad
I\_{\eta,i}=h_i^2\sum_jw\_{ij}s_j^2I\_{\mu,ij}.

The working target is \eta_i-g\_{\eta,i}/I\_{\eta,i}, with tiny
stabilization where the information vanishes. This quadratic
approximation updates the **log profile parameters**; the observed
intensities remain on their original scale and acceptance uses their
actual likelihood.

A smooth concave log profile can be written as

\eta(t)=a+bt-\frac12\sum_j c_j \\(t-u_j)\_+^2-(t-v_j)\_+^2\\,\qquad
c_j\geq0,

where the intervals \[u_j,v_j) partition the local time window. Its
second derivative is -\sum_jc_j\mathbf1\_{\[u_j,v_j)}(t)\leq0, giving
hard log-concavity. The intercept and slope are unrestricted, and a
small nonnegative quadratic solve estimates the curvature coefficients.
For more detailed local adjustments, piecewise-linear log-profile
increments preserve the existing fine-scale curvature; a
concave-envelope projection removes numerical constraint residuals.
Backtracking between feasible log profiles preserves concavity and
accepts only likelihood improvements.

The scaled active-set quadratic solver uses Cholesky factorizations for
positive-definite active subproblems. Extended-precision residuals
provide iterative corrections through the same factor. Conditioning and
residual checks trigger a rank-revealing QR fallback when needed. The
quadratic objective and shape constraints define the numerical
subproblem.

For a piecewise-linear log profile, each row of the interpolation matrix
has at most two nonzero entries. Its weighted Gram matrix is
tridiagonal, so natural-order sparse Cholesky and sparse
slope-constraint products avoid dense arithmetic on structural zeros.
Successive scoring steps initialize the inequality solver at its
preceding feasible solution and use a common KKT tolerance.

For a fixed component spectrum, a zero loading makes that channel’s
contribution constant in the profile update. Native workspaces therefore
retain its positive-loading channels in their original order. Spectral
updates and complete model comparisons still use every original channel.

Small joint conditional Poisson spectral and background refits use
observed-curvature Newton steps. Successive steps continue the
constrained active set, reuse fixed gradient scales and exploit sparse
weighted products for localized profiles. Conditional likelihood
backtracking verifies each step, with boundary checks and a KKT
tolerance controlling convergence.

The spectra act as covariates in each local profile update: Y\_{\mathcal
W}^{T}=S\_{\mathcal B}^{T}H\_{\mathcal W,\mathcal B}^{T} +R\_{\mathcal
W}^{T}. The current profiles and spectra provide the starts. Fine
refinement uses the constrained log-profile quadratic directly, with no
extra shrinkage toward zero or toward the initial PSF. Profile
normalization transfers its scale into the corresponding spectrum. Every
local block covers the full support of every component whose spectrum is
updated.

Fusion starts with spectral compatibility and profile overlap.
Broad-anchor and contained-tail proposals are refined only within the
affected windows; other compounds remain in an offset. The full
original-channel criterion decides whether each proposed rank change is
retained. A final conservative pass can reunite broad profiles with
complementary ion sets when their profile cosine is at least 0.95.
Neighbouring modes are guarded during its local refit.
`control=list(profile_cycles=0L)` disables the main shape-refinement
stage; `control=list(profile_fusion=FALSE)` disables the final
conservative pass.

A contained fragment can have a different profile from its broad parent
while carrying almost the same spectrum. The final contained-fragment
search therefore uses spectral cosine at least 0.97, at least 99% of the
fragment’s profile mass within the parent’s 1%-height interval, and a
fragment area no larger than 15% of the parent’s area. The parent FWHM
must exceed four times the representative PSF FWHM. A parent-profile
start and a log-concave combined-profile start are both tried. Each
candidate and an unmerged control receive six local joint-refinement
cycles, with four profile scoring steps per cycle. The complete
original-channel criterion chooses between them. This prevents
refinement gains alone from being attributed to a reduction in rank.

### Complete the fit on the original ion channels

The final pass checks conditional deletions, fusions of overlapping
profiles, and PSF candidates suggested by the original-channel
residuals. Nearly proportional profiles can be considered even when
neither is a broad peak. Candidate construction uses matrix products and
performs no refitting. Accepted proposals use native local likelihood
fits with all other compounds and the fitted background held in the
offset.

For J=-2\ell\_\nu+\lambda_RK, let \Delta
D=-2\\\ell\_\nu(M';Y)-\ell\_\nu(M;Y)\\. A deletion has change \Delta
J=\Delta D-\lambda_R, a fusion of r factors has \Delta J=\Delta
D-\lambda_R(r-1), and an addition has \Delta J=\Delta D+\lambda_R. The
complete original-channel criterion must improve. Noise parameters,
offsets, weights and the numerical compound penalty stay fixed. Cached
deletion scores need updating only for profiles whose temporal supports
overlap the deleted factor.

The pass uses bounded scan windows to screen missed peaks, then refines
the accepted profiles locally under hard log concavity. The same
procedure is available for every supported noise family; working-score
screening proposes candidates and each family’s likelihood decides
acceptance. The returned `fit$cleanup` records the moves and full
objective change. `control=list(original_cleanup=FALSE)` skips this
pass.

## Fusion, aggregation and rank selection

The estimated PSF supplies candidate shapes. Grouped sparse
deconvolution, fusion of spectrally similar overlapping components under
log-concavity, and local likelihood refinement jointly select a
parsimonious factorization rank. This yields a data-driven estimate of
the underlying number of components. Its distinctive feature is the
integration of these operations with an explicit noise-aware
compound-count objective.

1.  Obtain a variance-matched grouped deconvolution start, retaining the
    spline background as compulsory unpenalized terms.
2.  Compare normalized spectra by Jensen–Shannon information loss.
    Propose fusions when their elution profiles overlap and their
    combined profile can be made log-concave. Clustering itself performs
    no regression fits.
3.  Refine each proposal locally, holding all outside contributions
    fixed. Overlapping spectra and the candidate profile can adjust
    together.
4.  Accept a rank change when
    -2\\\ell\_\nu(M';Y)-\ell\_\nu(M;Y)\\+\lambda_R(K'-K)\leq0.
5.  Update broad profiles, smooth background coefficients and spectra,
    followed by final boundary checks and rank pruning.

Algebraic fusion uses a least concave majorant of the log profile,
followed by normalization. Local likelihood updates impose hard
concavity constraints in their quadratic solves. Every complete accepted
stage decreases the actual criterion, within numerical tolerance.
`history` and `ledger_error` make this property directly checkable.

With `precluster=32L`, the initial workflow uses at most 32 aggregated
outcome groups. Size factors `colSums(pmax(Y,0))` account for different
total ion intensities. \lambda_R retains the original channel count, and
final refinement and rank checks return to every original channel.
Poisson sums remain Poisson; common-dispersion quasipoisson sums retain
their mean-variance relationship. Independent Gaussian variances add.
Independent NB1 variables sharing \pi=1/\phi, with \phi\>1, are closed
under sums, because their probability-generating functions multiply as

\prod_j\left\\\frac{\pi}{1-(1-\pi)t}\right\\^{\mu_j/(\phi-1)}
=\left\\\frac{\pi}{1-(1-\pi)t}\right\\^{\sum_j\mu_j/(\phi-1)}.

The same property holds for common-scale Gamma latent signals, through
\prod_j(1+b_0t)^{-\mu_j/b_0}=(1+b_0t)^{-\sum_j\mu_j/b_0}. Standard NB2
and Gamma use moment matching in the reduced problem. For channel
intensity fractions r_j within a group G,

\theta_G=\theta/\sum\_{j\in G}r_j^2,\qquad \phi_G=\phi\sum\_{j\in
G}r_j^2,

respectively. These match the aggregate variance under the working
channel proportions. Censoring or different channel weights also give
working aggregate approximations; final scoring uses the original
responses and their specified model.

## Temporal variation

High-order differences suppress the slowly varying mean. With
a_k=(-1)^k{r\choose k} and stationary covariance c_h,

\frac{\operatorname{Var}(\Delta^r e)}{{2r\choose r}}
=c_0+2\sum\_{h=1}^r(-1)^h \frac{{2r\choose r+h}}{{2r\choose r}}c_h.

The identity follows by expanding the variance and using \sum_k a_k
a\_{k+h}=(-1)^h{2r\choose r+h}. With lag-one covariance this reduces to
c_0-2r c_1/(r+1). Fourth and eighth differences estimate the marginal
variance and lag-one covariance. The optional penalty multiplier
\kappa=1+2\widehat\rho_1 approximates the local information loss for
smooth directions, while the marginal scale stays in the component
likelihood. This is a scalar working covariance adjustment.
Composite-likelihood theory relates the broader correction to
\operatorname{tr}(JH^{-1}), where H is sensitivity and J is score
covariance (Varin, Reid & Firth, 2011).

## Single-sample algorithm

Let \mathcal L=-\ell\_\nu(M;Y). Resolve the pooled noise parameters \nu
and information multiplier \kappa before rank comparisons. The default
initializer uses a positive smoothed mean \widetilde\mu and
variance-matched Poisson working weights

w^{\rm pilot}\_{ij}=w\_{ij}\\
\frac{\widetilde\mu\_{ij}}{\operatorname{Var}(Y\_{ij}\mid\widetilde\mu\_{ij},\nu)}.

It fits grouped nonnegative `l0glmFit()` on the PSF design with
compulsory background, using the unadjusted \lambda\_{\rm IC}. The
starting search and original-channel polish use
J_0=-2\ell\_\nu+m\lambda\_{\rm IC}K. This permissive initialization
preserves candidates while the profiles are refined. Final rank
comparisons use J=-2\ell\_\nu+m\kappa\lambda\_{\rm IC}K. When \kappa=1,
the two criteria coincide. A supplied `start` with profiles, spectra and
optional background coefficients bypasses screening. The working Poisson
fit initializes the factorization; complete rank decisions use J with
the specified observation model.

**Algorithm N1 — sparse log-concave NMF with rank selection**

1.  Normalize candidate PSF columns. Resolve the noise scale, optional
    lag-one information adjustment, weights, offsets and background
    basis.
2.  If preclustering is requested:
    1.  Form size-adjusted outcome groups and sum their
        observations/offsets.
    2.  Construct the corresponding exact or moment-matched aggregate
        noise model.
    3.  Initialize and fuse components in that reduced space, retaining
        the original channel count in m\lambda\_{\rm IC}.
    4.  Lift loadings/background with group intensity fractions; return
        to every original channel for likelihood refinement.
3.  Otherwise obtain the variance-matched grouped deconvolution start
    directly. Normalize each profile and transfer its scale into its
    loadings.
4.  Run the bounded local refinement stages:
    1.  Hold other components and background fixed in a local offset.
        Update each spectrum by nonnegative scalar likelihood solves,
        including zero-boundary checks; update each concave log profile
        with the quadratic subroutine below.
    2.  Update nonnegative background coefficients conditionally on the
        component signal. Preserve the specified family likelihood.
    3.  Generate spectral/overlap-compatible fusion candidates
        algebraically, starting with stronger components. Refine only
        the affected scan windows.
    4.  Accept whole initialization candidates only if their complete
        J_0 decreases within numerical tolerance. Check conditional
        deletions and zero spectra.
    5.  If NB2 shape updating is enabled, propose bounded steps in
        \log\theta and backtrack until the original likelihood improves.
5.  On all original channels, run fine constrained profile/spectrum
    updates, conditional deletions, compatible fusions and
    residual-supported additions. Polish the accepted profiles again and
    revisit retained narrow-fragment candidates under J_0.
6.  Restore the requested information multiplier and evaluate the
    resulting initialization under J. Test contained fragments of broad
    components against equally refined unmerged local controls. Accept
    only complete original-channel improvements under this fixed final
    criterion.
7.  For nearly coincident modes, screen left-, right- and two-sided
    width restarts. Compatible Poisson/quasi-Poisson fits use the
    existing aggregated channels for this proposal screen; other cases
    use the original channels. Expand promising spectra using their
    current within-group proportions. Refine each candidate and an
    unperturbed control on all original channels, retaining a restart
    only when its complete J improves on both the current fit and the
    refined control. Reuse a control score for an identical local block
    only while the current fit is unchanged; invalidate it after a move.
8.  Return temporally ordered H,S, background, family/noise parameters,
    likelihood, rank criterion, stage acceptance ledger and timings.

`profile_cycles`, `fusion_passes` and the other stage controls bound the
work; local solvers also stop when their score/step change is
sufficiently small. The final restart pass is enabled by default when
`profile_cycles>0` and can be disabled with
`control=list(profile_restarts=FALSE)`.

For a profile h with mode a, define
T\_{l,r}(i)=a+\frac{i-a}{l\\\mathbf1\\i\<a\\+r\\\mathbf1\\i\geq a\\},
\qquad h_i'=\frac{\exp\\\eta(T\_{l,r}(i))\\}
{\sum_t\exp\\\eta(T\_{l,r}(t))\\},\quad\eta=\log h. The three starts use
(l,r)=(2,1),(1,2),(2,2), with piecewise-linear interpolation of \eta.
Concave slopes remain nonincreasing on each side, and their change from
nonnegative to nonpositive at the mode preserves log concavity across
the join. Candidate pairs have mode separation at most one quarter of
both FWHMs. Let P\in\\0,1\\^{m\times g} be the existing channel-group
assignment. For Poisson/quasi-Poisson with scalar dispersion and
likelihood weight and a zero offset, preliminary screening uses Y_a=YP,
S_a=SP and M_a=MP. Under the working independence model, summed counts
preserve the Poisson mean and the quasi-Poisson relationship
\operatorname{Var}(Y_a)=\phi M_a. No new clustering is needed.

For a promising candidate, expand each grouped spectral coefficient
using the entry model’s within-group ion proportions. Refine this
candidate and a separate unperturbed control on the **original
observations**. The final polish allows eight block cycles after
aggregated screening, or sixteen for the full-channel route, with local
convergence stopping earlier when appropriate. The control receives one
additional cycle. Noise, weights, offsets, background and component
count remain fixed, and acceptance requires

\ell\_\nu(Y;M\_{candidate})\>
\max\\\ell\_\nu(Y;M\_{entry}),\ell\_\nu(Y;M\_{control})\\.

Thus \Delta J=-2\Delta\ell\_\nu\<0 for an accepted restart. Aggregation
accelerates proposals; it does not replace the original-data objective
and need not choose the same local optimum as full-channel screening.
Other families, channel-specific scales, nonscalar observation weights
or nonzero offsets retain full-channel screening.

The cache stores only scalar control scores, not fitted mean matrices,
and is cleared after every accepted move. The `profile_starts` audit
records input component indices, starts, screening and original-data
changes, channel count, polish budget, cache hits/misses and time. Stage
timers measure elapsed work without forcing garbage collection between
nested fits; automatic R garbage collection remains active.

The initialization ledger records its unadjusted criterion. When the
final information multiplier differs from one, the final ledger begins
at the completed initialization and uses the fixed adjusted criterion
throughout. The noise scale remains on its marginal variance scale. The
reduced-space search supplies an initializer; final rank checks use
every original channel.

**Subroutine Q — local likelihood and shape update**

1.  Subtract the focus component from the current mean to obtain the
    fixed offset.
2.  Compute mean-scale likelihood scores and positive working
    information.
3.  For a loading b, propose a nonnegative scalar scoring step or use
    the bracketed Poisson root. For a profile, set \eta=\log h, apply
    the chain-rule scores and solve a weighted quadratic under
    nonincreasing slopes.
4.  Halve the step until the original likelihood improves and
    constraints hold. A convex combination of feasible log profiles
    remains concave.
5.  Normalize the profile, transfer its scale to loadings, and update
    only the affected mean block. Outside-window observations are
    unchanged.

For a rank change \Delta K, acceptance is based on

\Delta J=2(\mathcal L\_{\rm new}-\mathcal L\_{\rm old})+\lambda_R\Delta
K.

Because outside-window means are unchanged, their likelihood
contributions cancel exactly. Local computation therefore gives the
complete objective difference without rescoring unaffected observations.

## Simultaneous model on original scan grids

After suppressing backgrounds, exposures and offsets, the multisample
mean can be written M_s=Q_sP, where P is shared exactly and Q_s is
sample-specific. This is simultaneous NMF in the terminology of
[Cichocki et al.](https://doi.org/10.1002/9780470747278), [Badea
(2008)](https://psb.stanford.edu/psb-online/proceedings/psb08/badea.pdf)
and [Lee et al. (2012)](https://doi.org/10.1371/journal.pone.0048238).

For sample s, let u=w_s(t) be its monotone alignment map. Native bin
edges t\_{si-1/2},t\_{si+1/2} define

\Delta u\_{si}=w_s(t\_{si+1/2})-w_s(t\_{si-1/2}),\qquad
J\_{si}=\frac{\Delta u\_{si}}{\Delta t\_{si}}.

With known exposure e\_{si}\>0, the mean model is

\mu\_{sij}=e\_{si}\left\\J\_{si}\sum_g a\_{sg}h\_{sg}(u\_{si})p\_{gj}
+(B_sA_s)\_{ij}\right\\+o\_{sij}.

Spectra have p\_{gj}\geq0 and \sum_jp\_{gj}=1. Each detected profile
satisfies

h\_{sg}\geq0,\quad \sum_i\Delta u\_{si}h\_{sg}(u\_{si})=1,\quad
\frac{\log h\_{sg,i+1}-\log h\_{sg,i}}{u\_{s,i+1}-u\_{si}} \text{ is
nonincreasing}.

The Jacobian preserves integrated area under alignment. No interpolation
or resampling is applied to observations. Local likelihood evaluations
use their original channels, weights, offsets and sample-specific noise
parameters.

The rank-selection objective is

\mathcal J=-2\sum_s\ell\_{\nu_s}(M_s;Y_s) +\kappa\lambda\_{\rm
IC}\\(m-1)G+L+d\_\nu\\.

Here G is the number of common spectra and L the number of detected
sample-compound occurrences. A unit-sum spectrum contributes m-1 free
spectral coordinates and each occurrence contributes an amplitude.
d\_\nu=\sum_s d\_{\nu_s} counts estimated sample-specific noise
parameters, each shared across that sample’s channels. Profile and
baseline complexity are unpenalized shape-constrained nuisance terms in
this rank criterion. The weighted likelihood of each sample follows the
observation-model definition above, so observation weights are included
once.

For temporal-correlation adjustment, \kappa is the largest resolved
sample information factor. Quasipoisson uses each sample’s frozen,
scaled Poisson kernel. The same
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html),
[`EBIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
and
[`MBIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
interfaces apply.

For one sample, G=L=K, so

(m-1)G+L+d\_\nu=mK+d\_\nu.

At a common preset value and information factor, the simultaneous score
therefore reduces to the single-sample reported IC. Both select compound
count without an extra penalty on profile complexity or background
dimension.

### Penalty scale for single and simultaneous fits

| Interface | IC-name `lambda` | Numeric `lambda=c` |
|:---|:---|:---|
| [`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md) on one matrix | m\kappa\lambda\_{\rm IC} per compound | m\kappa c per compound |
| [`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md) on a sample list | \kappa\lambda\_{\rm IC} per selection degree of freedom | \kappa c per selection degree of freedom |

The simultaneous selection count is (m-1)G+L+d\_\nu. A compound observed
in one sample contributes m-1 spectral coordinates and one amplitude, so
its charge matches that of a single-sample component. The single-sample
solver
[`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
accepts a numeric total compound charge c; enter c/m in
[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
to reproduce that charge. Estimated noise parameters add a constant
selection count when their specification is fixed. Single-sample presets
use their original n,p; simultaneous presets use n\_\*=p\_\*=\sum_s n_s.

## Alignment and local updates

1.  Fit samples independently, optionally initializing in 32 aggregated
    channel groups and returning to every original channel.
2.  Match high-confidence spectra to a reference sample. Use the
    1%-height peak start plus half the narrow-peak width as the timing
    anchor.
3.  Fit a robust low-dimensional drift spline, with log area ratio as a
    nuisance term for concentration-related anchor motion. Retain its
    drift component and enforce positive slope.
4.  Group spectra algebraically when their timing and area-ordered
    profile overlap are compatible. Each group has at most one
    occurrence per sample. This clustering stage performs no NMF or GLM
    refits.
5.  Update local concave log profiles, pool spectral scores over
    samples, and update each smooth background. Backtrack against the
    original likelihood.
6.  Retain a grouped fit when its complete penalized objective improves.
    Return temporally ordered consensus spectra and sample-specific
    areas.

For a common spectral coefficient p\_{gj}, with design value
x\_{sgi}=e\_{si}J\_{si}a\_{sg}h\_{sg}(u\_{si}), the pooled
negative-log-likelihood score and positive working information are

U\_{gj}=\sum\_{si:g\in s}w\_{sij}x\_{sgi}g\_{\mu,sij},\qquad
I\_{gj}=\sum\_{si:g\in s}w\_{sij}x\_{sgi}^{2}I\_{\mu,sij}.

A nonnegative scoring step is followed by full-likelihood backtracking.
Renormalizing a spectrum transfers its scale into the occurrence
amplitudes and preserves the fitted signal. Profile updates use the
corresponding chain-rule score in log coordinates and linear slope
constraints on actual aligned times. Native OpenMP updates give
identical results across thread counts in the validation examples.

The native spectral kernel evaluates each cell’s likelihood derivatives
once and pools them over the sparse profiles in one pass parallelized
across channels. Cached peak geometry avoids rescanning every profile
after each algebraic merge. Shared row-wise weights and offsets retain
one vector throughout the simultaneous fit, with channel-specific
matrices used when required.

All seven [signal
families](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
are supported: Gaussian, Poisson, quasipoisson, NB2, NB1, censored Gamma
and common-scale Gamma. Dispersion parameters are resolved from the
original channels and frozen for group comparisons. NB2 additionally
supports likelihood-improving shape updates when initialized with a
negative `theta`.

## Simultaneous algorithm and normalization

For simultaneous models use n\_\*=p\_\*=\sum_s n_s in the IC preset, and
the largest sample information multiplier as \kappa. Noise parameters
are shared across channels within each sample. The joint rank count is
(m-1)G+L+d\_\nu, reducing to mK+d\_\nu when one sample has G=L=K. Common
unit-sum spectra contribute m-1 coordinates each; each observed
sample-compound occurrence contributes one amplitude.

If \sum_i\Delta u\_{si}h\_{sg}(u\_{si})=1, then

\sum_i\Delta t\_{si}J\_{si}h\_{sg}(u\_{si}) =\sum_i\Delta
u\_{si}h\_{sg}(u\_{si})=1.

Thus an occurrence amplitude represents integrated intensity in the
chosen time unit, before its known exposure is applied.
[`sampleNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
also returns unit-sum original-grid profiles and their compatible
discrete loadings.

**Algorithm N2 — simultaneous multisample log-concave NMF**

1.  Fit samples independently by N1, or read supplied fitted factors.
    Persist observations and factor checkpoints sample by sample.
2.  Match high-confidence spectra to a reference. Obtain
    onset-plus-half-width anchors and fit a robust monotone drift map,
    accounting for log area ratio as a nuisance predictor of anchor
    motion. Retain the drift component.
3.  Map each original bin’s edges to aligned time. Compute \Delta u,J,
    initialize log-concave profiles there, normalize by \Delta u, and
    put their scale into occurrence amplitudes. Initialize each
    occurrence with its own spectrum.
4.  Perform `initial_cycles` of local profile, pooled spectrum and
    background updates; retain this ungrouped-spectrum reference
    checkpoint.
5.  For up to `rounds` grouping rounds:
    1.  Restrict spectrum comparisons to nearby aligned time bins.
        Compare cosine or Hellinger distance, with at most one
        occurrence per sample per group.
    2.  Check area-ordered log-concave-ridge compatibility. Pool
        compatible spectra by occurrence area. Keep profiles fixed
        during this algebraic grouping.
    3.  Run up to `cycles` joint local-update cycles: constrained sample
        profiles; batched, pooled spectral scores; nonnegative sample
        backgrounds; optional NB2 shape updates. Backtrack each update
        on original-channel likelihood.
    4.  Stop grouping when a subsequent round makes no new joins.
6.  Compare the complete grouped/refined candidate against the reference
    checkpoint using the joint criterion. Retain the lower-criterion
    model.
7.  Compact consensus IDs in temporal order. Return shared spectra,
    sample factors, alignment, diagnostics and an area table with
    missing occurrences `NA`.

Clustering constructs a complete candidate before joint refinement. The
final comparison includes every original channel and all detected
occurrences.

### Area-ordered ridge compatibility

Put matched profiles on a small comparison grid for geometry only;
observed signals remain on their native grids. Order profiles by total
area and allow bounded local alignment. Their thresholded level sets
must overlap and support a common log-concave ridge across the area
ordering. This screens concentration- dependent broadening and tailing
before sharing a spectrum. The subsequent original-data likelihood
decides whether the combined proposal is retained.

### Streaming and parallel computation

Sparse profile supports restrict work to affected scan windows. The
simultaneous spectral kernel evaluates each observation’s likelihood
derivatives once and pools them across its active profiles. Channels can
be processed independently in OpenMP; sample contributions are
accumulated in a fixed order. The main optimization reads one sample at
a time, with shared spectra and compact geometry summaries in memory.
Shared row weights/offsets retain one vector. RT bins avoid an
all-occurrence pairwise similarity matrix.

## Methodological foundations and implementation

Sparse pulse-train regression builds on de Rooi & Eilers, and the
adaptive ridge initialization uses L0glm’s extensions of reweighted
quadratic selection. Eilers’ work motivates positive shape-constrained
profiles; generalized and deviance matrix factorization provide the
likelihood framework. Bregman and size-adjusted Poisson clustering
supply useful spectral and channel geometries.

The likelihood-based hard shape solves, complete original-channel rank
moves, windowed kernel checkpoints and simultaneous
normalization/acceptance rules above describe the logconNMF algorithms.
The simultaneous criterion counts shared spectra and sample-specific
amplitudes explicitly. The general references establish the underlying
frameworks; the objectives and pseudocode here specify their adaptation
to these constrained signal models.

## References

de Rooi, J. & Eilers, P. (2011). [Deconvolution of pulse trains with the
L0 penalty.](https://doi.org/10.1016/j.aca.2011.05.030)

Cichocki, A., Zdunek, R., Phan, A. H. & Amari, S.-I. (2009).
[*Nonnegative Matrix and Tensor Factorizations: Applications to
Exploratory Multi-way Data Analysis and Blind Source
Separation*.](https://doi.org/10.1002/9780470747278)

Badea, L. (2008). [Extracting gene expression profiles common to colon
and pancreatic adenocarcinoma using simultaneous nonnegative matrix
factorization.](https://psb.stanford.edu/psb-online/proceedings/psb08/badea.pdf)

Lee, C. M. et al. (2012). [Simultaneous non-negative matrix
factorization for multiple large scale gene expression datasets in
toxicology.](https://doi.org/10.1371/journal.pone.0048238)

Fan, Y. & Tang, C. Y. (2013). [Tuning parameter selection in high
dimensional penalized likelihood.](https://doi.org/10.1111/rssb.12001)

Ipsen, A. (2015). [Derivation from first principles of the statistical
distribution of the mass peak intensities of MS
data.](https://doi.org/10.1021/ac503554u)

Saumard, A. & Wellner, J. A. (2014). [Log-concavity and strong
log-concavity: a review.](https://doi.org/10.1214/14-SS107)

de Juan, A., Jaumot, J. & Tauler, R. (2014). [Multivariate curve
resolution (MCR). Solving the mixture analysis
problem.](https://doi.org/10.1039/C4AY00571F)

Kidziński, Ł., Hui, F. K. C., Warton, D. I. & Hastie, T. J. (2022).
[Generalized matrix factorization: efficient algorithms for fitting
generalized linear latent variable models to large data
arrays.](https://www.jmlr.org/papers/v23/20-1104.html)

Castiglione, C., Segers, A., Clement, L. & Risso, D. (2026). [Stochastic
gradient descent estimation of generalized matrix factorization models
with application to single-cell RNA sequencing
data.](https://doi.org/10.1093/biostatistics/kxag010)

Wang, L. & Carvalho, L. E. (2023). [Deviance matrix
factorization.](https://doi.org/10.1214/23-EJS2174)

Li, X. & Tyagi, A. (2023). [Block-active ADMM to minimize NMF with
Bregman divergences.](https://doi.org/10.3390/s23167229)

Banerjee, A. et al. (2005). [Clustering with Bregman
divergences.](https://jmlr.org/papers/v6/banerjee05b.html)

Witten, D. M. (2011). [Classification and clustering of sequencing data
using a Poisson model.](https://doi.org/10.1214/11-AOAS493)

Lee, D. D. & Seung, H. S. (2000). [Algorithms for non-negative matrix
factorization.](https://proceedings.neurips.cc/paper_files/paper/2000/file/f9d1152547c0bde01830b7e8bd60024c-Paper.pdf)

Eilers, P. H. C. (2005). [Unimodal
smoothing.](https://doi.org/10.1002/cem.935)

Varin, C., Reid, N. & Firth, D. (2011). [An overview of composite
likelihood
methods.](https://www3.stat.sinica.edu.tw/statistica/oldpdf/A21n11.pdf)

The likelihood distributions follow the R
[negative-binomial](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/NegBinomial.html)
and
[Gamma](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/GammaDist.html)
parameterizations. The compound-count criterion, local log-concave
fusion and their integration with L0glm’s grouped initialization are
developments of this package.

## Parallel disk-backed updates in the general identity-link solver

With common spectra P fixed, sample-specific factors are conditionally
independent. The simultaneous criterion is

J=-2\sum_s\ell_s(Y_s;\mu_s) +\lambda\\(m-1)G+L+d\_{\mathrm{noise}}\\,

where G is the number of shared compounds and L the number of detected
sample-compound occurrences. Profile and background updates leave this
rank charge fixed. They can therefore run concurrently over samples,
with each accepted update decreasing its own contribution to the full
criterion. Pooled spectral gradients and information are computed per
sample and summed in acquisition order. The shared spectral step is
accepted against the full original-channel likelihood.

A common channel map has a one-hot membership matrix
A\in\\0,1\\^{m\times q}, giving aggregated observations \bar Y_s=Y_sA
and spectra \bar P=PA. The map is learned once from summed time-bin
counts or supplied as a saved
[`channelGroups()`](https://tomwenseleers.github.io/logconNMF/reference/channelGroups.md)
object. Only this map-learning step uses a common time grid. The fitting
grids, bin exposures and sparse PSF matrices X_s remain sample specific.

Before matching, Poisson and common-scale quasipoisson observations with
shared scan weights use channel sums for profile proposals and all
original channels for acceptance. Once spectra are shared, every family
uses the native kernel of
[`fitLogConcaveProfile()`](https://tomwenseleers.github.io/logconNMF/reference/fitLogConcaveProfile.md)
on the original observations. For the focal component, the background
and overlapping neighbours form a fixed offset. A weighted projection of
the log rate onto the cone of nonincreasing slopes gives a scoring
direction on the aligned scan grid. A line search checks the resulting
mean against the family’s original-scale likelihood. The local window
extends two nominal narrow-peak FWHMs past each end of the initial
support and stays fixed across sweeps. Channels with zero focal loading
contribute a constant to this conditional objective and can be skipped.
Spectral, background and likelihood passes read at most 32 original
channels per block through a native binary-matrix reader. Thus
aggregation accelerates the initial proposal and matching steps while
the target criterion continues to use the original observations.

The execution order is:

1.  Learn or load the common channel map.
2.  Fit independent samples concurrently, each with its own sparse PSF.
3.  Estimate alignment and build the common spectral candidate graph.
4.  With spectra fixed, update profiles and backgrounds independently
    across samples; calculate pooled spectral contributions
    concurrently.
5.  Reduce contributions in sample order and update the shared spectra.
6.  Group compatible occurrences, refreshing only graph edges incident
    to a merged root. Distances and ordering of unaffected proposals are
    retained.
7.  Accept the complete candidate using the original-channel penalized
    objective and return full spectra, original-grid profiles and
    integrated areas.

`ncores` is a total logical-CPU budget; its automatic value for a sample
list is the detected logical core count. With W concurrent sample
workers, each receives \lfloor\texttt{ncores}/W\rfloor native threads.
BLAS runs with one thread during batch fitting and its caller setting is
restored afterwards. Only the assigned sample and compact factors are
sent to each worker; the complete observation cache is never serialized
to every worker.

During algebraic grouping, sparse profiles and aligned time vectors are
read once per sample into a transient geometry cache; this cache
contains no observations. Only merged roots need their occurrence groups
and area-weighted anchor centers updated. Rejected edges are skipped
when removed from the ordered queue, preserving the same greedy merge
order. Ridge integrals use the lower-area profile’s support and
bracketing zero grid points: all omitted trapezoidal contributions are
exactly zero. These changes reduce repeated file reads and work on
unrelated time points without changing spectral distances, compatibility
thresholds or original-channel acceptance criteria.

Eligible Poisson/quasipoisson workers retain at most one original matrix
of 64 MiB plus an equally sized matrix of \log\Gamma(y+1) normalizers.
The cache belongs to the active runtime, is replaced when the assigned
sample changes, and is released on success or error. Larger samples
retain local strided reads. Sixteen workers therefore add at most 2 GiB
of persistent original-data buffers, independent of the total batch
size. `control=list(cache=FALSE)` disables persistent input caches.
Reuse depends on sample affinity; a one-entry cache may be replaced
between passes when there are more samples than workers.

For a count observation the working negative log-likelihood is
L(y,\mu)=\\\mu-y\log\mu+\log\Gamma(y+1)\\/\phi, with \phi=1 for Poisson
noise. Gradients and Fisher information use (1-y/\mu)/\phi and
1/(\phi\mu), without evaluating the data-only normalizer. Absolute
likelihoods retain that constant through the cache; local likelihood
differences use the same stable difference kernel. A spectrum gradient
pass also supplies the baseline likelihood for its acceptance check.
Scaled spectra are reused within each local profile pass. These
execution details leave the original-channel objective and acceptance
thresholds intact.

Within coupled refinement, conditional Poisson backgrounds first check
the existing scaled KKT residual. Initialization and standalone
conditional refits retain their original Newton-update sequence. For
fixed peak means the background Hessian is
B^\mathsf{T}\operatorname{diag}(v_i y_i/\mu_i^2)B\succeq0, where v_i
includes observation weights and inverse dispersion. A stationary
nonnegative start therefore needs no further Newton solve. Nonstationary
starts retain the same constrained updates and convergence tolerance.
Noise-update passes are omitted when no sample has an iteratively
estimated NB2 shape.

The resident observation cache uses 8q\sum_s n_s bytes. For 1,000
samples with 10,000 scans and 32 groups, this is 2.56 GB, plus factors
and active-worker buffers. Original observations and factor checkpoints
are disk backed. Sparse sample-specific PSFs are modest in size: a
measured 14,946-by-14,946 design with 1,464,433 stored nonzero values
occupied 16.82 MiB. Sparse PSF workspace checks count unique structural
Gram entries and symbolic Cholesky fill; repeated overlaps are not
counted as distinct entries. Unpenalized spline columns are separated
before the Woodbury workspace check. The memory budget still applies to
designs with genuinely large Gram matrices or factor fill.
