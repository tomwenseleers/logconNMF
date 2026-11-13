# logconNMF 0.1.0.9000

* Corrected the historical one-unit mass-axis offset in `gcmsData1`: the
  all-zero leading channel was removed and the retained channels now have their
  nominal labels, m/z 30:699. `gcmsData2` was independently checked against its
  source CDF and required no correction.

## PSF estimation and deconvolution

* `convKernel()`, `shiftKernel()` and `convMatrix()` construct PSFs and
  convolution designs, including asymmetric bigaussian and EMG responses.
* `estimatePSF()` estimates constant or smoothly time-varying point spread
  functions from standards or overlapping signals and returns a sparse design.
* Alternating spike-train and kernel estimation supports single or multiple
  channels, smooth backgrounds and information-criterion checkpoints.

## Sparse log-concave matrix factorization

* `logconNMF()` estimates nonnegative profiles and loadings with hard
  log-concavity constraints and automatic best rank tuning.
* Local profile, spectrum and background updates, algebraic candidate fusion,
  conditional deletion and residual additions use original-channel criteria.
* Native profile fits use banded Cholesky systems with residual checks,
  iterative correction and a rank-revealing fallback. Cached local workspaces
  and continued active sets accelerate profile, spectral and background fits.
* Gaussian, Poisson, quasipoisson, NB2, NB1, censored Gamma and common-scale
  Gamma models, with pooled noise calibration and optional temporal adjustment.
* Optional channel aggregation accelerates initialization; final refinement
  retains every original outcome channel.
* Permissive initialization and native fine-scale profile/spectrum refinement
  preserve candidate components for final correlation-aware rank selection.
* Likelihood-checked local width restarts resolve signal allocation between
  nearly coincident profiles. Screening and acceptance retain all original
  channels and compare with an equally refined control.
* Contained broad-peak fragments are compared with equally refined unmerged
  controls. Complete local supports, reconstructed original-channel means
  and separate initialization/final-criterion ledgers make each decision auditable.

## Simultaneous multisample models

* `logconNMF()` with a sample list shares spectra across samples with sample-specific profiles,
  areas and backgrounds, allowing compounds to be absent from individual samples.
* Monotone alignment, area-ordered profile compatibility and exposure-aware
  likelihoods retain each sample's original irregular scan grid.
* Disk-backed observations, batched native spectral scores, shared row inputs,
  sparse profiles and OpenMP support efficient fitting.

* A total logical-core budget schedules independent fits and simultaneous sample
  updates concurrently, with deterministic spectral reductions and restored
  BLAS settings. Native block reads retain bounded original-channel buffers.
* Common saved channel maps, a compact aggregate cache and full-channel
  acceptance checks support reusable, memory-efficient multisample analyses.
* Incremental fusion candidates retain the same greedy order while refreshing
  only edges affected by a merge.

## Simulation and interpretation

* Single- and multisample simulations include drifting retention times,
  concentration-dependent width/tailing, broad profiles and missing compounds.
* Optional retention-index calibration places reference spectra along the
  chromatogram. Concentration-dependent broadening preserves leading edges,
  and each sample records its generating retention-time–RI map.
* Simulation baselines can be supplied as common levels, channel-specific
  intensities or smooth reference matrices on the original intensity scale.
* Experimental multichannel signals, profile and spectral plots, factor
  extraction, likelihood and information-criterion methods.
* Tutorials follow PSF estimation, single-sample factorization and simultaneous
  factorization; a separate mathematical article gives objectives and pseudocode.
