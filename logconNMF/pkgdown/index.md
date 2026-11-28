# logconNMF

<div class="l0-hero">
<div class="l0-hero-copy">
<p class="l0-eyebrow">SPARSE NONNEGATIVE MATRIX FACTORIZATION · LOG-CONCAVE PROFILES</p>
<h1>Log-concave NMF.<br><span>Automatic rank selection.</span></h1>
<p class="l0-lead">A structured sparse nonnegative factorisation framework for ordered multichannel data, combining automatic rank selection, hard log-concave profiles, flexible likelihoods supporting different noise models and simultaneous multisample models that do not require trilinearity.</p>
<div class="l0-actions"><a class="btn btn-primary" href="articles/log-concave-nmf.html">Get started <span aria-hidden="true">↗</span></a><a class="btn btn-outline-light" href="reference/index.html">Explore the methods</a></div>
<p class="l0-hero-note">R · Grouped L0 discovery · Penalised likelihood · Flexible multisample sharing</p>
</div>
<img class="l0-hero-logo" src="logo-dark.svg" alt="logconNMF hex logo with log-concave profiles and corresponding spectra" width="240" height="280">
</div>

<div class="l0-stat-grid">
<a href="reference/logconNMF.html"><strong>Automatic rank selection</strong><span>Score additions, deletions and fusions by penalised likelihood to obtain the most parsimonious rank</span></a>
<a href="articles/algorithm.html#log-concave-profiles"><strong>Hard log-concavity</strong><span>Component- and sample-specific chromatographic profiles</span></a>
<a href="articles/log-concave-nmf.html"><strong>Flexible noise families</strong><span>Gaussian, Poisson, quasi-Poisson, Negative Binomial and Gamma models</span></a>
</div>

## Nonnegative matrix and tensor factorisation to separate overlapping peaks

Hyphenated chromatography is a demanding test case for structured NMF and
tensor decomposition. In GC–MS, chromatographic separation is rarely
complete: when compounds co-elute, the spectrum observed at one scan is an
additive mixture of their spectra. In matrix notation,

$$
\mathbf{Y}_{i\boldsymbol{\cdot}}
\approx \sum_{r=1}^{R} H_{ir}\mathbf{S}_{r\boldsymbol{\cdot}}
+ \mathbf{B}_{i\boldsymbol{\cdot}},
$$

where $\mathbf{H}_{\boldsymbol{\cdot}r}$ is the elution-mode factor,
$\mathbf{S}_{r\boldsymbol{\cdot}}$ is the spectral-mode factor for component
$r$, and $\mathbf{B}_{i\boldsymbol{\cdot}}$ represents the smooth background
at scan $i$. Recovering the latent components is therefore a
constrained nonnegative factorisation problem: model order must be identified
while chemically meaningful structure is imposed on both modes.

Nonnegative matrix and tensor factorisation have a substantial history in
hyphenated chromatography. Important examples include
[MCR-ALS](https://doi.org/10.1039/C4AY00571F),
[NMF applied directly to GC–MS](https://doi.org/10.1016/j.talanta.2004.09.017),
[nonnegative unimodal matrix factorisation](https://doi.org/10.1109/ICASSP39728.2021.9414631),
[global parametric analysis of multiple GC/MS data sets](https://doi.org/10.1016/j.chemolab.2008.10.004),
[PARAFAC2 for retention-time shifts](https://doi.org/10.1002/%28SICI%291099-128X%28199905/08%2913%3A3/4%3C295%3A%3AAID-CEM547%3E3.0.CO%3B2-Y),
and the PARAFAC2-based
[PARADISe workflow](https://doi.org/10.1016/j.chroma.2017.04.052), including its recent
[shift-invariant trilinearity formulation](https://doi.org/10.1002/cem.3501).
These methods established that multichannel structure can resolve co-elution
and recover spectra with less interference than apex extraction alone.
For the factorisation community, these data expose where unconstrained
bilinear models and rigid multilinear coupling become fragile: rank, profile
geometry, observation noise, background and between-sample variation all
interact.

### What standard formulations leave unresolved

- **Optimisation and chemical rank.** Bilinear and multilinear factorisations
  are nonconvex. Different starts can give different stationary solutions,
  while ALS algorithms may enter long *swamps* or stop at a poor local fit
  ([Boutsidis & Gallopoulos, 2008](https://doi.org/10.1016/j.patcog.2007.09.010);
  [Tomasi & Bro, 2006](https://doi.org/10.1016/j.csda.2004.11.013)). The number
  of components is usually supplied or selected outside the main fit. Useful
  diagnostics and automated classifiers exist for PARAFAC/PARAFAC2
  ([Bro & Kiers, 2003](https://doi.org/10.1002/cem.801);
  [Johnsen et al., 2014](https://doi.org/10.1002/cem.2575)), but require fitting
  and comparing a sequence of prespecified ranks.
- **Peak geometry.** Nonnegativity alone permits multimodal or jagged elution
  profiles. MCR-ALS can impose unimodality
  ([Bro & Sidiropoulos, 1998](https://doi.org/10.1002/%28SICI%291099-128X%28199807/08%2912%3A4%3C223%3A%3AAID-CEM511%3E3.0.CO%3B2-2)),
  but unimodality by itself permits shoulders, plateaus and irregular tails.
  A more appropriate peak-shape constraint would be log-concavity, which is stronger:
  when imposed on the positive support,

  $$
  \Delta^2\log h_i
  = \log h_{i+1}-2\log h_i+\log h_{i-1}\leq 0,
  $$

  which implies unimodality while still allowing asymmetric and long-tailed
  peaks.
- **Instrument response, background and noise.** Standard factorisations do
  not jointly estimate a constant or smoothly time-varying point-spread
  function (PSF) and a smooth nonnegative background. Least squares also
  imposes a constant-variance Gaussian working model. Ion production and
  arrival motivate count models, while detector gain, aggregation and chemical
  variability can produce overdispersion or approximately Gaussian intensities
  ([Ipsen, 2015](https://doi.org/10.1021/ac503554u)). Although Poisson MCR has
  been proposed—in the NMF literature, the corresponding Poisson-likelihood
  formulation is known as Kullback–Leibler NMF
  ([Cichocki et al., 2009](https://doi.org/10.1002/9780470747278);
  [Lavoie et al., 2016](https://doi.org/10.1016/j.chemolab.2016.02.006))—it is a
  specialised extension that has not been routinely used in the chemometrics
  literature and does not account for overdispersion.
- **Multiple samples.** PARAFAC assumes trilinearity: apart from an amplitude,
  each compound has the same elution profile in every sample. Shifted PARAFAC
  adds translations but not compound- and concentration-specific broadening or
  tailing. Classical PARAFAC2 is more flexible, but still requires
  $\mathbf{H}_s^{\mathsf T}\mathbf{H}_s=\boldsymbol{\Phi}$ for every sample
  $s$. Retention shifts can preserve
  these cross-products; changing widths, tails and pairwise overlaps generally
  do not. Flexible and fully constrained PARAFAC2 variants relax important
  parts of this restriction
  ([Cohen & Bro, 2018](https://doi.org/10.1007/978-3-319-93764-9_9);
  [Roald et al., 2022](https://doi.org/10.1137/21M1450033)).

::: {.l0-method-comparison}

| Method (typical cited formulation) | Rank selection | Chromatographic profile | Nonnegativity | Joint component selection across channels | Optimal spectral sparsity | Noise objective | Multiple samples | Rank and scope in the cited formulation |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| [MCR-ALS / NMF](https://doi.org/10.6028/jres.124.018)<sup>†</sup> ([GC–MS NMF](https://doi.org/10.1016/j.talanta.2004.09.017); [unimodal NMF](https://doi.org/10.1109/ICASSP39728.2021.9414631)) | Supplied; external EFA/PCA diagnostics | Optional unimodality and other constraints, but no hard log-concavity | Both factors in standard NMF and NNLS-constrained MCR; optional in broader MCR | No | No | Least squares is usual; generalized KL and specialised Poisson extensions exist | Separate, augmented or shared-factor multiset fits | Chosen for the mixture or interval being resolved |
| PARAFAC / shifted PARAFAC | Supplied; external diagnostics | One shared profile per compound; shifted variant adds translation | Optional | No | No | Least squares in standard algorithms | Exact trilinearity, or sample shifts only | Chosen for the sample tensor or interval |
| [Classical PARAFAC2 / PARADISe](https://doi.org/10.1016/j.chroma.2017.04.052) | Usually supplied; multi-rank diagnostic automation has been proposed | Sample-specific profiles with a constant cross-product constraint; no hard chromatographic shape | Historically not on every mode | No | No | Least squares | Retention drift and restricted shape variation | PARADISe compares candidate ranks within each chromatographic interval |
| Flexible / [fully constrained PARAFAC2](https://doi.org/10.1137/21M1450033) | Supplied externally | Optional nonnegativity, smoothness or unimodality; relaxed/flexible coupling | Can constrain every mode | No | No | Least squares in the cited formulation | More flexible drift and shape variation, with an explicit coupling penalty or constraint | Fitted at the chosen order |
| [Shift-invariant trilinearity](https://doi.org/10.1002/cem.3501) | Supplied externally | A shared peak shape up to translation; cannot model sample-specific shape changes | Can constrain all modes | No | No | Least squares | Retention shifts under shift-invariant trilinearity | Evaluated on low-rank windows; whole chromatograms still require many separate models |

:::

<p class="l0-table-footnote"><sup>†</sup> When both factors are constrained nonnegative and each alternating least-squares step uses nonnegative least squares under a squared-error loss, MCR-ALS and standard least-squares NMF solve the same constrained bilinear factorisation (<a href="https://doi.org/10.6028/jres.124.018">Camp, 2019</a>). <strong>MCR-ALS</strong> remains the established name in chemometrics and is the broader framework because its regressors and constraints need not be nonnegative least squares.</p>

The cited chromatographic workflows choose rank for a matrix or retention interval. In PARADISe, intervals are selected to contain preferably fewer than six peaks. The rank of a whole chromatogram counts compounds across all intervals. [PARADISe workflow](https://doi.org/10.1016/j.chroma.2017.04.052)

### Commercial GC–MS software

AnalyzerPro exposes peak-width and detection thresholds. MassHunter Unknowns Analysis can search multiple widths, with a trade-off between resolving close peaks and runtime. LECO ChromaTOF Sync co-processes Pegasus sample sets to improve spectra and align peaks. GCMSprime fits one shared spectrum and a separate log-concave profile per occurrence, and scores candidate splits against the original ion counts. The commercial workflows have yet to be scored on these 16 chromatograms. [AnalyzerPro](https://spectralworks.com/2016/07/08/using-analyzerpro-in-the-metabolomics-workflow/) · [MassHunter](https://www.agilent.com/cs/library/usermanuals/public/G3335-90243_Unknowns_Analysis_Familiarization-en.pdf#page=27) · [ChromaTOF Sync](https://www.leco.com/instrument/chromatof-sync-data-alignment/) [Nine commercial workflows compared](https://tomwenseleers.github.io/GCMSprime/articles/benchmarks.html#commercial-gcms-deconvolution-software).

### What logconNMF adds to constrained factorisation

<div class="l0-feature-grid">
<div><span class="l0-feature-number">01 · GROUPED DISCOVERY</span><h4>Grouped peak selection</h4><p>Candidate peak positions are selected jointly across detector channels, pooling evidence across mass channels rather than thresholding ions separately.</p></div>
<div><span class="l0-feature-number">02 · SPECTRAL SPARSITY</span><h4>Optimal spectral sparsity</h4><p>Within retained components, sparsity-promoting penalised regression favours parsimonious ion support; final selection moves are checked on every original channel.</p></div>
<div><span class="l0-feature-number">03 · MODEL SIZE</span><h4>Automatic rank selection</h4><p>Score additions, deletions and fusions by penalised likelihood to obtain the most parsimonious rank.</p></div>
<div><span class="l0-feature-number">04 · PEAK GEOMETRY</span><h4>Hard log-concave profiles</h4><p>Every retained profile is nonnegative, contiguous and log-concave during fitting—a stronger constraint than post-hoc unimodality.</p></div>
<div><span class="l0-feature-number">05 · INSTRUMENT RESPONSE</span><h4>PSF starts, individual final shapes</h4><p><code>estimatePSF()</code> estimates a constant or time-varying PSF for initialisation. The constant blind-PSF model is a constrained lag-separable convolutive NMF; each final compound occurrence nevertheless receives its own width and tail.</p></div>
<div><span class="l0-feature-number">06 · OBSERVATION MODEL</span><h4>Background and noise fitted jointly</h4><p>A smooth nonnegative background remains in the model. Flexible likelihood families cover count-like, overdispersed and continuous measurements.</p></div>
<div><span class="l0-feature-number">07 · MULTIPLE SAMPLES</span><h4>Shared spectra without trilinearity</h4><p>Samples retain their original grids, shapes, locations and areas. Monotone area-preserving time warps account for drift without resampling, and compounds may be absent.</p></div>
<div><span class="l0-feature-number">08 · WHOLE-DATA SCALE</span><h4>High-rank full chromatograms</h4><p>Sparse screening and local likelihood updates support full chromatograms and large batches directly, with hundreds or thousands of components or sample–component occurrences.</p></div>
</div>

Structured PSF starts, grouped screening, original-channel refinement, profile
restarts and criterion-checked fusion reduce sensitivity to poor starts. The
time-varying PSF is locally rather than globally convolutional, and the final
factorisation does not constrain profiles to shifted copies of a common kernel.
The same formulation applies to other ordered multichannel signals with
nonnegative, single-peaked components.

### Rank-selection objective

For ordinary identity-link count data, the original-ion objective is

$$
J = \sum_s D_{\rm Pois}(Y_s,\mu_s)/\phi_s +
\lambda_{\rm eff}\left[\sum_k(\|p_k\|_0-1)+\sum_{s,k}d(h_{sk})\right].
$$

Each shared unit-sum spectrum contributes one fewer degree of freedom than
its number of positive ions. Each sample-specific log profile contributes the
dimension of its fitted piecewise-linear shape. Absent profiles contribute
zero. The same criterion checks changes in peak shape, ion support and
component count for one or several samples. The fixed background basis stays
unpenalized. Coordinate updates and local support exchanges can reach a local
minimum. The [mathematical description](articles/algorithm.html#conditional-count-objective) gives the acceptance rules
and the criterion used for other noise families.

## Installation

```r
install.packages("remotes")
remotes::install_github(c(
  "tomwenseleers/L0glm/L0glm",
  "tomwenseleers/logconNMF/logconNMF"
))
```

The R package includes the conditional Poisson and quasi-Poisson solver. It fits the original ion counts with compiled C++ kernels and uses OpenMP when available. Solver selection is automatic for the requested noise model and memory settings.

## Log-concave NMF for single samples

This deliberately small cuticular-hydrocarbon window contains six co-eluting
components and is used to make the recovered factors easy to inspect. It is
not a scale limit: logconNMF runs efficiently on full sample matrices with
approximately 10,000–15,000 time scans or more. The initial dictionary `X` places a narrow
Gaussian point spread function (PSF) at each possible peak location. An
efficient [L0glm](https://tomwenseleers.github.io/L0glm/) fit then selects the
best subset that can explain the multichannel signal, thereby penalising most
coefficients to zero. Equally, logconNMF provides
[`estimatePSFFromStd()`](reference/estimatePSFFromStd.html) and
[`estimatePSF()`](reference/estimatePSF.html) to estimate the PSF and construct
the most parsimonious corresponding dictionary from either externally run
standards, such as an alkane ladder, or the sample chromatogram itself.

```r
library(logconNMF)
data(gcmsData1)
Y <- gcmsData1
g <- convKernel("gaussian", 65L, list(sigma = 9))
psf <- estimatePSF(Y, time = as.numeric(rownames(Y)), standard = NULL,
  nwindows = 1L, precluster = 1L,
  control = list(kernel_est = g, maxiter = 0L))
X <- psf$X
fit <- logconNMF(Y, X = X, time = as.numeric(rownames(Y)), background = matrix(1, nrow(Y), 1L),
  family = quasipoisson("identity"), precluster = 7L,
  control = list(initial_relaxed = FALSE, fusion_passes = 0L))

fit$rank                         # estimated number of components
profiles <- fit$H                 # scans x components
loadings <- coef(fit)             # components x original channels
spectra <- coef(fit, normalize = "l1")
GIC(fit)
```

Using `estimatePSF()` here applies the same boundary-aware, column-normalised
PSF construction as the high-level GCMSprime workflow. Consequently, these
settings reproduce the identical fitted profiles, spectra and background in
both packages; repeated fits are deterministic.

The columns of `profiles` sum to one, so each row sum in `loadings`
gives the component's total fitted signal across scans and channels.
Grouping similar channels can speed up the initial search. The reported
spectra use every original channel. We selected this six-component example
after a sensitivity grid over Gaussian PSF widths and precluster counts: among
solutions with six factors and chemically plausible retention times, it gave
the best balance of low GIC, accurate apex times, compact profiles and strong
enrichment of the compounds' diagnostic ions.

![Background-corrected ion traces in grey, six fitted log-concave profiles and a 2 by 3 panel of recovered spectra for six co-eluting components](reference/figures/gcms-data1-decomposition.png)

The upper square-root display makes weak ion traces visible; the six spectra
below are shown as base-peak relative intensities.

### Beer volatiles: compound-specific log-concave peak shapes

A beer-headspace window provides a more demanding real-data example. Narrow,
partially overlapping volatiles occur alongside a broad, strongly asymmetric
acid peak. Instead of forcing one common PSF on every compound, the final
factors act as compound-specific log-concave effective peak-shape functions:
they can be narrow, broad or strongly asymmetric while remaining nonnegative
and single-peaked. The associated spectra are estimated jointly from all 518
mass channels.

[![Beer-headspace chromatogram with fitted compound-specific log-concave profiles](reference/figures/gcms-data2-traces.png)](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html#beer-volatiles-and-a-broad-acid-peak)

Observed ion traces are shown in grey; coloured curves show the fitted
compound contributions. Component 17 captures the broad asymmetric acid signal
without imposing its shape on the narrower neighbouring peaks.

## Simultaneous log-concave NMF for multiple samples

For related matrices, logconNMF performs **simultaneous NMF** in the terminology
of [Cichocki et al. (2009)](https://doi.org/10.1002/9780470747278),
[Badea (2008)](https://psb.stanford.edu/psb-online/proceedings/psb08/badea.pdf)
and [Lee et al. (2012)](https://doi.org/10.1371/journal.pone.0048238):

$$
\mathbf{Y}_s
\approx \mathbf{H}_s\mathrm{diag}(\mathbf{a}_s)\mathbf{S}
+\mathbf{B}_s\mathbf{A}_s.
$$

Here $\mathbf{H}_s$ contains unit-area log-concave profiles,
$\mathbf{a}_s$ gives the component areas, $\mathbf{S}$ contains the shared
unit-sum component spectra, and $\mathbf{B}_s\mathbf{A}_s$ is the smooth
nonnegative background. Monotone area-preserving time warps model retention
drift and are taken into account as alignment and Jacobian offsets without
resampling samples on a common grid; components may be absent from any subset
of samples. This shares chemical
information without imposing PARAFAC trilinearity or the PARAFAC2
constant-cross-product constraint.

Grouped PSF regressions divide the core budget between independent samples and native threads. The scheduler estimates how many regressions fit in available memory, including Windows commit headroom. After those workers exit, a single process shares the original ion counts and fitted means. OpenMP updates elution profiles across samples. It fits each ion’s shared spectrum and sample-specific background coefficients independently, keeping a fixed component order within each ion. The general identity-link solver can stream original-channel blocks from disk.
### Automatic rank recovery on full chromatograms

The same 16 simulated chromatograms were supplied to every program. Each contains **10,000 scans × 550 nominal mass channels** (m/z 30–579), covering approximately 50 minutes. Together they contain **88 million simulated ion intensities**. Each of the **1,000 distinct compounds** occurs in eight randomly chosen samples. The design therefore contains **8,000 positive sample–compound pairs**. Rank here counts the distinct generating compounds across the batch, each with a time-local elution profile in each positive sample. Of the compounds, 990 have exponentially modified Gaussian (EMG) peaks and ten have a rising shark-fin profile with an abrupt trailing edge. Higher concentrations broaden both types. The acquisition clock allows a smooth retention shift of up to 2.93 seconds between the first and last injections and approximately 2% variation in scan spacing.

GCMSprime, using simultaneous logconNMF, **recovered the true rank of 1,000 exactly** in **35.95 seconds per sample**. It recovered **8,000 of 8,000** compound occurrences and made **no false-positive calls**. The recovered peaks included **80 of the 80 acid occurrences**. MSHub returned **1,414 deconvolved components**, 41.4% above the true rank. The other five workflows produce peak lists or aligned feature groups. Their outputs do not define a common factorisation rank.

The 1,000 shared factors matched the 1,000 generating compounds one-to-one across all eight occurrences of each compound. This is exact automatic recovery of the chemical component count from 550 measured ion channels. Scan order distinguishes compounds that elute at different times, and each compound contributes a shared spectrum wherever it occurs.

Published rank-selection experiments put the scale in perspective. Distributed NMF (DnMFk) recovered seven planted factors in a one-terabyte matrix. Eswar et al. compared NMF rank selectors at true rank 10. Fu et al. tested nonnegative tensor rank selection at true ranks 15–40. The 1,000 compounds here are **25 times the largest generating rank in that tensor experiment**. These studies use different data models and assess different kinds of rank. [Chennupati et al.](https://doi.org/10.1007/s11227-020-03181-6) · [Eswar et al.](https://www.osti.gov/servlets/purl/3002147) · [Fu et al.](https://doi.org/10.3390/e21070632)

| Workflow | Seconds/sample | Estimated rank (true 1,000) | TP | FP | FN | Acid TP / 80 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GCMSprime | 35.95 | **1,000 (exact)** | 8,000 | 0 | 0 | 80 / 80 |
| MS-DIAL | 16.02 | — | 7,728 | 11,685 | 272 | 73 / 80 |
| mzmine GC-EI | 6.92 | — | 4,202 | 185 | 3,798 | 29 / 80 |
| eRah | 25.94 | — | 7,172 | 5,615 | 828 | 70 / 80 |
| AMDIS | 28.99 | — | 7,878 | 9,122 | 122 | 80 / 80 |
| xcms + CAMERA | 4.22 | — | 7,261 | 890 | 739 | 31 / 80 |
| MSHub | 15.05 | 1,414 | 7,179 | 12,184 | 821 | 50 / 80 |

Among GCMSprime’s 8,000 matched occurrences, the median absolute relative area error was **0.095%** and ordinary adjusted R² on the log₁₀ scale was **0.999982**.

[![Estimated and true areas in the simultaneous fit](reference/figures/simultaneous-benchmark-area-accuracy.png)](reference/figures/simultaneous-benchmark-area-accuracy.png)

[![Accuracy of recovered shared spectra](reference/figures/simultaneous-benchmark-spectral-accuracy.png)](reference/figures/simultaneous-benchmark-spectral-accuracy.png)

[![Recovery of ions within each spectrum](reference/figures/simultaneous-benchmark-spectral-sparsity.png)](reference/figures/simultaneous-benchmark-spectral-sparsity.png)

### Optimality gap to ground truth-initialised solution

The blind fit showed smooth, monotonic descent in GIC across eight sweeps. Its final GIC differed by only **0.19%** from that of a fit initialised with the true elution profiles and run for the same eight sweeps. Both fits ended at the true rank of **1,000** with **8,000** sample-level occurrences.

[![Original-ion GIC over eight sweeps of the blind fit on 16 full chromatograms](reference/figures/full-benchmark-gic-gap.png)](reference/figures/full-benchmark-gic-gap.png)

[Block-level GIC trace](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/gic-trace.csv) · [Final GIC values](https://github.com/tomwenseleers/GCMSprime/blob/master/GCMSbenchmark/inst/benchmarks/results/gic-gap.csv).

### Elution profile recovery

The examples are the EMG compound and the acid with the greatest total true area in their respective classes. Both were selected from the simulation before comparing any fits. Each occurs in eight samples. Solid curves show inferred profiles, dashed curves show the true profiles, and the eight absent true profiles coincide at zero. Colours order the 16 samples by true area. Heights retain original total-ion intensity per scan. The mean cosine in each example uses all eight positive occurrences, with a missing recovered profile scoring zero. The all-acid tally covers all ten acids.

[![EMG and shark-fin acid profiles in all 16 samples](reference/figures/simultaneous-benchmark-elution-profiles.png)](reference/figures/simultaneous-benchmark-elution-profiles.png)

[Benchmark methods, scoring, timings and reproduction code](https://tomwenseleers.github.io/GCMSprime/articles/benchmarks.html).

## Choose a workflow

| Data and objective | Method |
|:--|:--|
| Estimate a common response shape from a standard or overlapping signals | [Estimate PSFs](articles/estimate-psf.html) with `estimatePSF()`. The returned sparse design can use a constant or time-varying PSF. |
| Resolve components with individual shapes in one matrix | [Sparse log-concave NMF](articles/log-concave-nmf.html) with `logconNMF()`. The tutorial covers rank selection, backgrounds and flexible noise families. |
| Share spectra across related samples | [Simultaneous log-concave NMF](articles/simultaneous-nmf.html) with `logconNMF()`. Each sample retains its own scan grid, peak shapes and areas. Components may be absent from individual samples. |

For a known and common PSF, [L0glm](https://tomwenseleers.github.io/L0glm/articles/deconvolution.html)
fits the sparse convolution model directly. [GCMSprime](https://tomwenseleers.github.io/GCMSprime/)
uses logconNMF for chromatogram files, retention-index calibration and library
search. The [algorithm and theory](articles/algorithm.html) give the
objectives, shape constraints and mathematical pseudocode.

## Selected literature

The comparison above draws on the following methodological and
chromatography-specific work.

**Matrix factorisation and curve resolution**

- de Juan A, Jaumot J, Tauler R (2014).
  [Multivariate Curve Resolution (MCR). Solving the mixture analysis problem.](https://doi.org/10.1039/C4AY00571F)
- de Juan A, Tauler R (2021).
  [Multivariate Curve Resolution: 50 years addressing the mixture analysis problem—A review.](https://doi.org/10.1016/j.aca.2020.10.051)
- Camp CH Jr (2019).
  [pyMCR: A Python Library for Multivariate Curve Resolution Analysis with Alternating Regression (MCR-AR).](https://doi.org/10.6028/jres.124.018)
- Gao H-T, Li T-H, Chen K, Li W-G, Bi X (2005).
  [Overlapping spectra resolution using non-negative matrix factorization.](https://doi.org/10.1016/j.talanta.2004.09.017)
- Ang AMS, Gillis N, Vandaele A, De Sterck H (2021).
  [Nonnegative Unimodal Matrix Factorization.](https://doi.org/10.1109/ICASSP39728.2021.9414631)
- van Stokkum IHM, Mullen KM, Mihaleva VV (2009).
  [Global analysis of multiple gas chromatography–mass spectrometry (GC/MS) data sets: A method for resolution of co-eluting components with comparison to MCR-ALS.](https://doi.org/10.1016/j.chemolab.2008.10.004)
- Lavoie FB, Braidy N, Gosselin R (2016).
  [Including noise characteristics in MCR to improve mapping and component extraction from spectral images.](https://doi.org/10.1016/j.chemolab.2016.02.006)

**Simultaneous and convolutive NMF**

- Cichocki A, Zdunek R, Phan AH, Amari S-I (2009).
  [*Nonnegative Matrix and Tensor Factorizations: Applications to Exploratory Multi-way Data Analysis and Blind Source Separation*.](https://doi.org/10.1002/9780470747278)
- Badea L (2008).
  [Extracting gene expression profiles common to colon and pancreatic adenocarcinoma using simultaneous nonnegative matrix factorization.](https://psb.stanford.edu/psb-online/proceedings/psb08/badea.pdf)
- Lee CM, Mudaliar MAV, Haggart DR, Wolf CR, Miele G, Vass JK, Higham DJ (2012).
  [Simultaneous non-negative matrix factorization for multiple large scale gene expression datasets in toxicology.](https://doi.org/10.1371/journal.pone.0048238)

**Tensor factorisation for chromatographic data**

- Bro R, Andersson CA, Kiers HAL (1999).
  [PARAFAC2—Part II. Modeling chromatographic data with retention time shifts.](https://doi.org/10.1002/%28SICI%291099-128X%28199905/08%2913%3A3/4%3C295%3A%3AAID-CEM547%3E3.0.CO%3B2-Y)
- Amigo JM, Skov T, Coello J, Maspoch S, Bro R (2008).
  [Solving GC–MS problems with PARAFAC2.](https://doi.org/10.1016/j.trac.2008.05.011)
- Johnsen LG, Amigo JM, Skov T, Bro R (2014).
  [Automated resolution of overlapping peaks in chromatographic data.](https://doi.org/10.1002/cem.2575)
- Johnsen LG, Skou PB, Khakimov B, Bro R (2017).
  [Gas chromatography–mass spectrometry data processing made easy.](https://doi.org/10.1016/j.chroma.2017.04.052)
- Schneide P-A, Bro R, Gallagher NB (2023).
  [Shift-invariant tri-linearity—A new model for resolving untargeted gas chromatography coupled mass spectrometry data.](https://doi.org/10.1002/cem.3501)
- Cohen JE, Bro R (2018).
  [Nonnegative PARAFAC2: a flexible coupling approach.](https://doi.org/10.1007/978-3-319-93764-9_9)
- Roald M, Schenker C, Calhoun VD, Adalı T, Bro R, Cohen JE, Acar E (2022).
  [An AO-ADMM approach to constraining PARAFAC2 on all modes.](https://doi.org/10.1137/21M1450033)
- Sorochan Armstrong MD, Hinrich JL, de la Mata AP, Harynuk JJ (2023).
  [PARAFAC2×N: coupled decomposition of multi-modal data with drift in N modes.](https://doi.org/10.1016/j.aca.2023.340909)

**Shape constraints, rank diagnostics and statistical modelling**

- Bro R, Sidiropoulos ND (1998).
  [Least squares algorithms under unimodality and non-negativity constraints.](https://doi.org/10.1002/%28SICI%291099-128X%28199807/08%2912%3A4%3C223%3A%3AAID-CEM511%3E3.0.CO%3B2-2)
- Eilers PHC (2005).
  [Unimodal smoothing.](https://doi.org/10.1002/cem.935)
- Bro R, Kiers HAL (2003).
  [A new efficient method for determining the number of components in PARAFAC models.](https://doi.org/10.1002/cem.801)
- Boutsidis C, Gallopoulos E (2008).
  [SVD based initialization: a head start for nonnegative matrix factorization.](https://doi.org/10.1016/j.patcog.2007.09.010)
- Tomasi G, Bro R (2006).
  [A comparison of algorithms for fitting the PARAFAC model.](https://doi.org/10.1016/j.csda.2004.11.013)
- de Rooi JJ, Eilers PHC (2011).
  [Deconvolution of pulse trains with the L0 penalty.](https://doi.org/10.1016/j.aca.2011.05.030)
- Ipsen A (2015).
  [Derivation from first principles of the statistical distribution of the mass peak intensities of MS data.](https://doi.org/10.1021/ac503554u)

## Citation

An article describing the methods implemented in logconNMF is in preparation.
If you use the package before the article is formally published, we would
appreciate a citation to the package:

Wenseleers T, Vanderaa C (2026).
*logconNMF: Log-Concave Nonnegative Matrix Factorization with Automatic Rank Selection*.
R package version 0.1.0.9000,
<https://tomwenseleers.github.io/logconNMF/>.

The [Authors and Citation page](authors.html#citation) also provides the
citation in BibTeX format. Please also cite the methodological sources linked
in the help pages where relevant.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize and
refactor code and help draft the documentation. The key ideas and methodological
decisions are the authors' own, and the authors take full responsibility for the
work presented here.
