# Sparse deconvolution with a fixed PSF via L0-penalised GLMs

Estimates the *latent*, un-blurred single or multichannel spike train
\\a\\ from an observed signal \\y\\ that is assumed to follow the
convolution model \$\$y \\=\\ X(g)\\a + \varepsilon,\$\$ where \\X(g)\\
is a Toeplitz design matrix that contains time-shifted copies of a known
point-spread function \\g\\ (PSF, a.k.a. convolution kernel) and errors
\\\varepsilon\\ are assumed to be drawn from an exponential family error
distribution (e.g. Poisson, Gaussian, Gamma, Inverse Gaussian). The
coefficient vector (or matrix) \\a\\ is found with *best-subset* or
*best-subset of group selection* via L0glm.

## Usage

``` r
l0glmSparseDeco(
  y,
  family = poisson(identity),
  kernel = NULL,
  kernelfunc = function(x) convKernel("gaussian", L = length(x), params = list(sigma =
    6), perc = 0),
  L = if (!is.null(kernel)) length(kernel) else 51,
  X = NULL,
  baseline = NULL,
  nknots = max(4L, floor(nrow(as.matrix(y))/500)),
  IC = "gic",
  grouped = TRUE,
  control.solver = do.call("controlSolverGen", list(solver = "eigen", method = 1L)),
  offset = NULL,
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

- kernel:

  Either a numeric vector that is the fixed PSF in forward time order,
  or \`NULL\` together with \`L\` and \`kernelfunc\`, in which case the
  PSF is generated on the fly (useful for simulations).

- kernelfunc:

  A function that returns a length-\\L\\ numeric vector when called as
  \`kernelfunc(1:L)\`. Defaults to a Gaussian with sigma=6.

- L:

  Desired kernel length when \`kernel\` is \`NULL\`.

- X:

  Optional. A sparse design matrix that overrides the automatic call to
  [`convMatrix`](https://tomwenseleers.github.io/logconNMF/reference/convMatrix.md).
  This is useful when you need time-varying PSFs.

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

- IC:

  Information criterion for the penalty preset; see
  [information_criteria](https://tomwenseleers.github.io/L0glm/reference/information_criteria.html).

- grouped:

  Logical, passed unchanged to `l0glmFit`, specifying whether to perform
  grouped L0 feature selection, i.e. best subset of group selection
  (BSGS), in case there are multiple outcome channels. Defaults to TRUE.

- control.solver:

  Optional solver-control list passed to `l0glmFit`. By default this
  uses Eigen sparse Cholesky LLT with AMD ordering. For strictly banded
  PSF designs, users can pass
  `controlSolverGen("eigen", method=1L, NaturalOrdering=TRUE)` to
  request Eigen natural ordering.

- offset:

  Optional known additive signal: a scalar, length-\\n\\ vector, or
  \\n\times m\\ matrix. Original observations are fitted with this model
  offset; counts are not background-subtracted.

- ...:

  Execution controls passed to the spike/background
  [`l0glmFit()`](https://tomwenseleers.github.io/L0glm/reference/l0glmFit.html)
  calls, such as `ncores` and `nthreads`. Built-in scalar families use
  native fitting automatically, with `ncores` limiting OpenMP
  outcome-channel threads. Set `ncores=1L` for serial fitting. Custom
  callbacks remain serial and use this budget for BLAS unless `nthreads`
  overrides it.

## Value

A list with

- `coef_spikes`:

  Estimated un-blurred spikes, i.e. the coefficients for design matrix
  `X_spikes`. Its dimensions are `ncol(X_spikes)` by
  `ncol(as.matrix(y))`.

- `coef_baseline`:

  Estimated baseline coefficients for baseline design matrix
  `X_baseline`.

- `fitted_values`:

  \\X a\\, the fitted signal.

- `fitted_baseline`:

  \\X a\\, the fitted baseline

- `deviance`:

  Residual deviance of the fit.

- `X_spikes`:

  Toeplitz design matrix for the spikes.

- `X_baseline`:

  Design matrix for the baseline.

- `offset`:

  Known additive signal, expanded to the response dimensions.

- `fit`:

  The underlying `L0glm_fit` model, supporting
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`GIC()`](https://tomwenseleers.github.io/L0glm/reference/model_criteria.html)
  and the other fitted-model methods.

## Details

The algorithm is inspired by the L0-penalised deconvolution approach of
de Rooi & Eilers (2011) for chromatographic pulse trains, but it
leverages the general machinery of L0glm. This results in the following
improvements:

- No manual tuning of the regularisation parameter is required; instead
  it can be preset to approximate information criteria such as GIC, AIC,
  BIC, EBIC or MBIC

- Both single channel data with L0 best subset selection or multichannel
  data with grouped L0 best subset of group selection are allowed for

- Any exponential-family noise model with an *identity* link can be used
  (Poisson, Gaussian, Gamma, ...)

- The underlying sparse design is handled by Eigen's efficient sparse
  Cholesky solver, enabling large problems to be tackled quickly

In the non-negative-matrix-factorisation (NMF) literature this model
would be referred to as *single-atom convolutional NMF* or *sparse
convolutional dictionary learning*, whereby a one-dimensional signal (or
multichannel matrix) \\y\\ is approximated as \\W(g)\\H\\, where
\\W(g)\\ is a banded Toeplitz matrix that holds time-shifted copies of a
fixed dictionary atom/PSF \\g\\, and \\H\\ is a (potentially
multi-column) matrix of sparse, non-negative activations \\a\\. However,
relative to standard single-atom convolutional NMF we offer several
extensions:

- *Optimal sparsity*: we do not fix the sparsity level of the
  activations \\a\\ in advance; instead we let L0glm perform best-subset
  or best subset of group selection through an L\\\_{0}\\ or grouped
  L\\\_{0}\\ penalty,

- *General noise models*: any exponential-family likelihood with an
  *identity* link is supported (Poisson, Gaussian, Gamma,
  inverse-Gaussian, ...), instead of only Gaussian (Euclidean-distance
  NMF) or Poisson (NMF with KL divergence)

- *Speed*: the Toeplitz design is solved by sparse-Cholesky linear
  algebra rather than multiplicative updates, giving orders-of-magnitude
  speed-ups on long signals

Methodologically this could be referred to as *generalised / deviance
matrix factorisation* (Kidzinski et al. 2022; Castiglione et al. 2024;
Wang & Carvalho 2023), specialised to the case where \\W\\ is a Toeplitz
single-atom convolutional dictionary matrix and where the matrix
factorisation is solved with an L\\\_{0}\\ or grouped
L\\\_{0}\\-penalised generalised linear model (GLM) fit via adaptive
ridge regularised iteratively reweighted least squares (IRLS). When a
fitted baseline or known offset is present, positive backprojection
starts are calculated from the signal above a low channel-wise baseline.
This keeps promising peaks eligible during the initial adaptive-ridge
steps. These values initialize the coefficients only; the model always
fits the original observations with its fitted background and known
offset.

A classical baseline for Poisson deconvolution is the Richardson-Lucy
algorithm (Richardson 1972; Lucy 1974). In modern terms, Richardson-Lucy
equals Kullback-Leibler NMF with multiplicative updates (Lee & Seung
1999) and factor matrix \\W\\ constrained to be a convolution matrix.
Empirically, `l0glmSparseDeco()` has much better support recovery,
attains much lower reconstruction error and is dramatically faster
thanks to (i) closed-form sparse linear solves, (ii) the automatic
selection of the penalty parameter through information-criterion
presets, and (iii) the possibility to exploit grouped-sparsity when
several channels are deconvolved jointly.

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

[`l0glmDecoEstKernel`](https://tomwenseleers.github.io/logconNMF/reference/l0glmDecoEstKernel.md)
for the complementary "estimate-the-kernel" step,
[`l0glmBlindDeco`](https://tomwenseleers.github.io/logconNMF/reference/l0glmBlindDeco.md)
for the full alternating blind deconvolution,
[`convMatrix`](https://tomwenseleers.github.io/logconNMF/reference/convMatrix.md)
(design matrix) and
[`convKernel`](https://tomwenseleers.github.io/logconNMF/reference/convKernel.md)
(PSF generator).

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
n        <- 500L                    # timepoints
npeaks   <- 50L                     # nonzero peaks
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
## 2  Deconvolution with the known PSF
## ------------------------------------------------------------------------

fit <- l0glmSparseDeco(y, kernel = g_true,
                         family = poisson(identity),
                         baseline = baseline, 
                         IC = "gic")
a_est_baseline <- as.vector(fit$coef_baseline) 
a_est <- as.vector(fit$coef_spikes)
sum(a_est!=0) 
#> [1] 52

## ---------- summary statistics ------------------------------------
n_true <- sum(a_true != 0)
n_est  <- sum(a_est  != 0)
tp     <- sum(a_true != 0 & a_est != 0)
precision <- tp / n_est
recall    <- tp / n_true
F1        <- 2 * precision * recall / (precision + recall)
R2p <- 100 *
  summary(lm(a_true ~ 0 + a_est))$r.squared
R2p
#> [1] 98.23545

## ---------- scatter-plot with true vs. estimated peak areas -------
main_expr <- bquote(
  atop("Poisson L0-glm deconvolution with known PSF",
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


# Fit several outcome channels jointly, with their own baseline coefficients.
Y <- cbind(y, rpois(n, 0.7 * y_true), rpois(n, 1.3 * y_true))
multi <- l0glmSparseDeco(Y, kernel = g_true, family = poisson(identity),
                          baseline = "bspline", nknots = 6L)
dim(multi$coef_spikes)
#> [1] 500   3
dim(multi$coef_baseline)
#> [1] 6 3
```
