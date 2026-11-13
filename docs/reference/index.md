# Package index

## Estimate point spread functions and deconvolve signals

- [`convKernel()`](https://tomwenseleers.github.io/logconNMF/reference/convKernel.md)
  [`shiftKernel()`](https://tomwenseleers.github.io/logconNMF/reference/convKernel.md)
  : Generate normalised 1-D convolution kernels (point-spread functions)
- [`convMatrix()`](https://tomwenseleers.github.io/logconNMF/reference/convMatrix.md)
  : Make banded or rectangular convolution design matrices
- [`estimatePSF()`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSF.md)
  : Estimate constant or time-varying point spread functions
- [`estimatePSFFromStd()`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSFFromStd.md)
  : Estimate a time-varying PSF matrix from an alkane ladder TIC
- [`l0glmBlindDeco()`](https://tomwenseleers.github.io/logconNMF/reference/l0glmBlindDeco.md)
  : Alternating sparse L0 glm powered blind deconvolution
- [`l0glmSparseDeco()`](https://tomwenseleers.github.io/logconNMF/reference/l0glmSparseDeco.md)
  : Sparse deconvolution with a fixed PSF via L0-penalised GLMs
- [`l0glmDecoEstKernel()`](https://tomwenseleers.github.io/logconNMF/reference/l0glmDecoEstKernel.md)
  : Estimate the PSF given a spike train estimate via L0-penalised GLMs

## Sparse and simultaneous factorization

- [`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
  [`coef(`*`<logconNMF_fit>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
  [`fitted(`*`<logconNMF_fit>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
  [`plot(`*`<simultaneousNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
  [`sampleNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
  : Fit a log-concave factorization with automatic rank selection
- [`residuals(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`splogconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`fitted(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`nobs(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`coef(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`logLik(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`AIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`BIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`GIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`EBIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`MBIC(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  [`print(`*`<splogconNMF>`*`)`](https://tomwenseleers.github.io/logconNMF/reference/splogconNMF.md)
  : Sparse log-concave nonnegative matrix factorization
- [`channelGroups()`](https://tomwenseleers.github.io/logconNMF/reference/channelGroups.md)
  [`aggregateChannels()`](https://tomwenseleers.github.io/logconNMF/reference/channelGroups.md)
  : Learn and reuse a common channel aggregation

## Noise calibration and profile constraints

- [`estimateSignalNoise()`](https://tomwenseleers.github.io/logconNMF/reference/estimateSignalNoise.md)
  : Estimate a shared noise model for ordered multichannel signals
- [`logConcSmooth()`](https://tomwenseleers.github.io/logconNMF/reference/logConcSmooth.md)
  : Log-concave spline smoothing of a discrete convolution kernel
  estimate

## Simulations and experimental data

- [`simulateGCMS()`](https://tomwenseleers.github.io/logconNMF/reference/simulateGCMS.md)
  : Simulate synthetic GC-MS chromatograms with ground-truth signal +
  background
- [`simulateGCMSSamples()`](https://tomwenseleers.github.io/logconNMF/reference/simulateGCMSSamples.md)
  : Simulate related chromatograms with shared compound spectra
- [`gcmsData1`](https://tomwenseleers.github.io/logconNMF/reference/gcmsData1.md)
  : GC/MS Chromatogram Section - Cuticular Hydrocarbon Profile of a
  Cataglyphis iberica Ant Queen
- [`gcmsData2`](https://tomwenseleers.github.io/logconNMF/reference/gcmsData2.md)
  : GC/MS Chromatogram Section - Beer Headspace Volatiles

## Shared tools and package overview

- [`reexports`](https://tomwenseleers.github.io/logconNMF/reference/reexports.md)
  [`GIC`](https://tomwenseleers.github.io/logconNMF/reference/reexports.md)
  [`EBIC`](https://tomwenseleers.github.io/logconNMF/reference/reexports.md)
  [`MBIC`](https://tomwenseleers.github.io/logconNMF/reference/reexports.md)
  : Objects exported from other packages
- [`logconNMF-package`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF-package.md)
  : Log-Concave Nonnegative Matrix Factorization with Automatic Rank
  Selection
