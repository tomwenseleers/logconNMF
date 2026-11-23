# logconNMF 0.1.0.9008

- `truth.profiles` compares a conditional count fit with a second fit initialized from matched simulation profiles. `gicGap()` reports the signed empirical GIC gap and its plot shows both trajectories.

# logconNMF 0.1.0.9007

- The conditional count solver runs within R using compiled C++ and OpenMP. It uses the same penalized original-ion objective for profile, spectral and local support updates. The general identity-link solver supports the other noise families, weights, offsets and exposures.
- Sparse retention profiles provide the interpolated 1%-height onsets for alignment. Reciprocal spectral comparisons use the retention window directly. Original-ion spectral regression, local offsets and support changes are evaluated in compiled code.

# logconNMF 0.1.0.9006

- Multisample count fits can resolve a mixed spectrum into several independently observed spectra. Overlapping profiles are refitted on original channels, with each change checked against the same penalized objective.

# logconNMF 0.1.0.9005

* The joint count solver uses the requested core budget for independent sample profiles and background updates, and for ion columns of each shared spectrum. Spectral coordinates retain their update order within each ion. Original-channel GIC and log-concavity checks also apply to threaded fits.

# logconNMF 0.1.0.9004

* The default Poisson and quasi-Poisson count fit alternates original-ion updates of hard log-concave profiles and sparse spectra. The GIC counts fitted profile and spectrum dimensions. A single matrix and registered sample lists use the same criterion. `objective.trace=TRUE` records the accepted objective values. `algorithm="classic"` selects the general-family implementation.
* Grouped PSF regressions in large batches cap their worker count using available memory. The conditional solver releases the worker pool before refining the factors.

# logconNMF 0.1.0.9003

* `splogconNMF()` and single-sample `logconNMF()` accept `control$initial_solver`,
  a named list of `L0glm::controlSolverGen()` arguments for component selection.
* Cholesky workspace-budget failures during initialization automatically retry
  with diagonally preconditioned Eigen CG and report the switch. The retry
  preserves the workspace budget, data and selection settings. Unrelated errors
  propagate; `control$initial_fallback=FALSE` disables the retry.
* `cache=FALSE` also disables temporary whole-sample observation caches in
  coupled profile refinement, retaining local-window disk reads.

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
  BLAS settings. Automatic sample-list budgets use logical cores, while single
  samples use physical cores. Native block reads retain bounded original-channel buffers.
* Common saved channel maps, a compact aggregate cache and full-channel
  acceptance checks support reusable, memory-efficient multisample analyses.
* Incremental fusion candidates retain the same greedy order while refreshing
  only edges affected by a merge.
* Coupled refinement reuses baseline likelihoods, local scaled spectra and
  derivative-only kernels. Eligible count workers retain one original matrix
  up to 64 MiB plus its likelihood normalizers; larger samples stream blocks.
  Runtime-scoped caches are released on success or error, and `cache=FALSE`
  disables persistent input buffers. Absolute scores keep all likelihood
  constants, and original-channel acceptance checks are retained.
* Coupled background updates recognize starts satisfying the existing KKT
  tolerance; initialization retains its original Newton-update sequence;
  unchanged noise parameters avoid redundant passes. Iteration budgets,
  grouping rules and the statistical objective are unchanged.
* Algebraic grouping caches sparse profiles without observations, updates only
  changed root centers, and restricts ridge integrals to their nonzero local
  supports while preserving complete-grid scores and greedy decisions.

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
