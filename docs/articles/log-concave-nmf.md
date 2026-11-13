# Sparse log-concave NMF with automatic best rank tuning

A GC–MS chromatogram can contain several compounds at the same retention
time. The time course of a compound is its elution profile, and the
relative signals it produces across ions form its mass spectrum. Each
measured scan contains a mixture of those contributions.

A fixed convolution model places the same local response shape at each
candidate position. That works well when the peaks have similar shapes.
Broad, tailed or overloaded peaks can require several adjacent spikes to
approximate one compound.
[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
begins with such a dictionary and then estimates a separate log-concave
profile and spectrum for every retained component. This lets the model
combine compatible fragments and refine both the peak shapes and the
spectra.

## The factorization

For n ordered scans and m channels, write the mean as

M = H S + B A + O,\qquad H\in\mathbb{R}\_+^{n\times K},\quad
S\in\mathbb{R}\_+^{K\times m}.

Here K is the number of components, BA is an optional smooth background
and O is a known additive offset. Component k contributes the outer
product of column k of H and row k of S. The positive support of a
profile is contiguous. On an evenly spaced grid, log-concavity means
2\log H\_{ik}\geq \log H\_{i-1,k}+\log H\_{i+1,k} within that support.
This includes familiar Gaussian and exponentially tailed peak shapes.

Each column of H sums to one. The corresponding row of S then contains
integrated channel intensities and can be normalized to a spectrum of
relative intensities. The fit compares likelihood and rank so that
another profile must explain enough of the signal to justify its
inclusion.

## Fit profiles, loadings and background

The signal below contains six co-eluting components. Rows index ordered
scans; columns index mass channels. The formulation also applies to
other ordered multichannel measurements.

``` r

library(L0glm)
library(logconNMF)
data(gcmsData1)
Y <- gcmsData1
X <- convMatrix(convKernel("gaussian", 65L, list(sigma=8)), nrow(Y))
B <- matrix(1, nrow(Y), 1L)
fit <- logconNMF(Y, X=X, background=B,
  family=quasipoisson("identity"), precluster=8L, correlation="auto", ncores=8L)
fit$rank
dim(fit$H)
dim(coef(fit))
GIC(fit)
```

The constant background is supplied as a single column of ones. Use [PSF
estimation](https://tomwenseleers.github.io/logconNMF/articles/estimate-psf.md)
to build `X` when the initial response shape is unknown or varies with
position.

## Inspect the recovered components

``` r

matplot(as.matrix(fit$H), type="l", lty=1, xlab="Scan", ylab="Profile")
matplot(t(coef(fit)), type="h", lty=1, xlab="Channel", ylab="Loading")
L0glm::plotMat(t(fitted(fit)), xlab="Scan", ylab="Channel")
patterns <- coef(fit) / rowSums(coef(fit))
fitted_signal <- fitted(fit)
fit$seconds
```

![Background-corrected ion traces, six resolved log-concave profiles and
a 2 by 3 panel of recovered
spectra](../reference/figures/gcms-data1-decomposition.png)

Background-corrected ion traces, six resolved log-concave profiles and a
2 by 3 panel of recovered spectra

All 670 channels are retained, and factors are returned in increasing
profile-apex order. The upper panel uses a square-root display to retain
weak grey ion traces; the six spectra use base-peak relative
intensities. The [GCMSprime
tutorial](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
develops the chromatographic interpretation, file input and reports.

## A broad asymmetric component

The beer-volatile window combines narrow overlapping peaks with a much
broader acid contribution. A low-dimensional spline accounts for the
smooth background:

``` r

data(gcmsData2)
Y2 <- gcmsData2
t2 <- as.numeric(rownames(Y2))
X2 <- convMatrix(convKernel("gaussian", 65L, list(sigma=6)), nrow(Y2))
B2 <- splines::bs(t2, df=4L, intercept=TRUE)
fit2 <- logconNMF(Y2, X=X2, background=B2,
  family=quasipoisson("identity"), precluster=32L, correlation="auto", ncores=8L)
fit2$rank
fit2$contained_fusion$broad
```

![Log-concave profiles and full-channel spectra, including the broad
component](../reference/figures/gcms-data2-profiles-spectra.png)

Log-concave profiles and full-channel spectra, including the broad
component

The table reports one direct fit of each matrix with eight OpenMP
threads and one BLAS thread on an Intel Core i9-11950H. Constructing the
PSF dictionary and background basis is outside the timing.

| Window             | Scans | Channels | Components | Quasi-GIC | Seconds |
|:-------------------|------:|---------:|-----------:|----------:|--------:|
| Hydrocarbon window |   119 |      670 |          6 |   247,888 |    7.96 |
| Beer volatiles     | 1,764 |      518 |         38 | 1,467,426 |   67.00 |

The [GCMSprime
tutorial](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)
times the complete CDF-to-MSP workflow separately and provides
diagnostic plots and downloadable spectra.

## Automatic best rank tuning

The initialization keeps a permissive pool of components while native
constrained likelihood updates refine their profiles and spectra. Final
rank selection applies the requested temporal information adjustment.
Small fragments contained within a broad profile are compared with an
equally refined unmerged local model on all original channels.

The default `lambda="gic"` compares complete fits by likelihood and
component count. After dividing the likelihood by the original channel
count, each retained component incurs one rank charge. On the unscaled
likelihood the charge is m\lambda\_{\mathrm{IC}}, with an optional
temporal information adjustment. Profile complexity and background
coefficients have no extra rank penalty.

``` r

fit$objective
fit$history[, c("stage", "delta", "accepted", "compounds")]
fit$cleanup
fit$initialization
fit$contained_fusion
c(AIC=AIC(fit), BIC=BIC(fit), GIC=GIC(fit))
```

Reported criteria count mK, where K is the fitted rank, plus estimated
shared noise parameters. Other presets are `"aic"`, `"bic"`, `"ebic"`
and `"mbic"`. A nonnegative numeric `lambda` supplies the penalty per
selection degree of freedom. Each retained component incurs m\lambda
before temporal adjustment. The [rank
objective](https://tomwenseleers.github.io/logconNMF/articles/algorithm.html#rank-selection)
gives the full convention.

## Choose a noise model

All seven models use identity-link nonnegative means and the original
observation scale. Gaussian measurements can have negative errors around
that mean.

| Model | R specification | Variance at mean \mu |
|:---|:---|:---|
| Gaussian | `gaussian("identity")` | \phi |
| Poisson | `poisson("identity")` | \mu |
| Quasipoisson | `quasipoisson("identity")` | \phi\mu |
| NB2 | `mgcv::nb(link="identity")` | \mu+\mu^2/\theta |
| NB1 | `nb1(dispersion=2, link="identity")` | \phi\mu |
| Censored Gamma | `gammaCensored(dispersion=.2, link="identity")` | \phi\mu^2 before censoring |
| Common-scale Gamma | `gammaScale(scale=2, link="identity")` | \sigma\mu before censoring |

Here \phi denotes the dispersion and \sigma the common Gamma scale.
Positive `dispersion` or `theta` arguments fix the parameter. Omitted
shared parameters are estimated from temporal differences. A negative
NB2 `theta` initializes at its absolute value and enables
likelihood-improving updates. Quasipoisson uses a scaled Poisson kernel
and corresponding quasi-criteria.

Zero Gamma intensities represent left censoring at a detection limit,
set by `control=list(censor=0.5)`; their contribution is the probability
below that limit. See the [noise
models](https://tomwenseleers.github.io/logconNMF/articles/algorithm.html#noise-models)
for scores and parameterizations.

## Faster initialization and smooth backgrounds

`precluster=32L` initializes in at most 32 size-adjusted outcome groups,
then returns to every original channel. `precluster=FALSE` uses the full
matrix directly. `control$linkage` selects Ward Poisson or Bregman
grouping. The original channel count sets the rank charge in both
workflows; final scoring uses the original observations and specified
noise model.

For a longer signal, replace the constant background with a
low-dimensional nonnegative spline basis:

``` r

B <- splines::bs(seq_len(nrow(Y)), df=4L, intercept=TRUE)
```

Background coefficients are fitted conditionally on the component signal
and remain unpenalized in rank selection. Known additive signals use
`offset`; weights can be shared across channels or channel-specific.

The default `correlation="auto"` estimates a lag-one working information
adjustment. `"none"` uses independent-observation rank charges; a
supplied nonnegative correlation fixes the adjustment. Marginal noise
variance stays in the likelihood. The [mathematical
article](https://tomwenseleers.github.io/logconNMF/articles/algorithm.md)
develops local solvers, hard shape constraints, aggregation and complete
pseudocode. Continue with [simultaneous multisample
NMF](https://tomwenseleers.github.io/logconNMF/articles/simultaneous-nmf.md)
for shared loadings and sample-specific profiles.
