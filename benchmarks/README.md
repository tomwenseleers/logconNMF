# Reproduce the examples

Run `Rscript benchmarks/benchmark_gcms_windows.R` with the installed package
on the R library path to fit both bundled matrices directly. The
[single-sample tutorial](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.html)
reports elapsed fitting time, component counts and quasi-GIC. The
[simultaneous-NMF tutorial](https://tomwenseleers.github.io/logconNMF/articles/simultaneous-nmf.html)
gives the three-sample simulation and links to the complete CDF benchmark.

The direct fits used R 4.5.0 on an Intel Core i9-11950H with one BLAS thread,
logconNMF 0.1.0.9011 and L0glm with response-scale IRLS initialization.
The scripts record thread budgets and package versions alongside the results.
