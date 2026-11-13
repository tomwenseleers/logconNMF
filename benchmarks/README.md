# Reproduce the examples

Run `Rscript benchmarks/benchmark_gcms_windows.R` with the installed package
on the R library path to fit both bundled matrices directly. The
[single-sample tutorial](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.html)
reports elapsed fitting time, component counts and quasi-GIC. The
[simultaneous-NMF tutorial](https://tomwenseleers.github.io/logconNMF/articles/simultaneous-nmf.html)
gives the three-sample simulation and links to the complete CDF benchmark.

The November 2026 runs used logconNMF
[`1c85840`](https://github.com/tomwenseleers/logconNMF/commit/1c85840)
and L0glm [`d7bdece`](https://github.com/tomwenseleers/L0glm/commit/d7bdece),
with R 4.5.0 and one BLAS thread on an Intel Core i9-11950H workstation.
The tables give the OpenMP thread or worker budgets for each run.
