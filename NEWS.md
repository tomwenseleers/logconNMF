# logconNMF 0.1.0.9004

* The default Poisson and quasi-Poisson count fit alternates original-ion updates of hard log-concave profiles and sparse spectra. The GIC counts fitted profile and spectrum dimensions. A single matrix and registered sample lists use the same criterion. `objective.trace=TRUE` records the accepted objective values. `algorithm="classic"` selects the general-family implementation. See the mathematical article for installation of the optional Python numerical backend.
* Grouped PSF regressions in large batches cap their worker count using available memory. The conditional solver releases the worker pool before refining the factors.

# logconNMF 0.1.0.9003

* Expose component-initialization solver controls via `control$initial_solver`.
  Cholesky workspace-budget errors retry with diagonally preconditioned Eigen
  CG, preserving the memory limit and selection settings and reporting the
  switch. Disable the retry with `control$initial_fallback=FALSE`.
* Disabling coupled observation caches also prevents temporary full-sample
  profile caches, retaining local-window disk reads.

# logconNMF 0.1.0.9002

* Bundled `gcmsData1` channels use nominal m/z 30:699; experimental mass axes
  are checked against the source data.

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
  nearly coincident profiles. Compatible Poisson/quasi-Poisson proposal screens
  reuse channel aggregates; completed candidates and matched controls always
  use every original channel. Scalar control-score caches are invalidated after
  each accepted move.
* Contained broad-peak fragments are compared with equally refined unmerged
  controls. Complete local supports, reconstructed original-channel means
  and separate initialization/final-criterion ledgers make each decision auditable.
* Nested stage timing avoids forced garbage collection between local fits;
  ordinary R memory management remains active.

## Simultaneous multisample models

* `fitLogConcaveProfile()` estimates a rate profile for a fixed spectrum,
  background and exposure under any of the seven identity-link noise models.
  After cross-sample matching, coupled fits use this update on the original
  channels and evaluate shared spectra with the joint information criterion.

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

* Coupled refinement reuses likelihoods and bounded worker caches, with
  unchanged objectives, grouping rules and original-channel checks.
