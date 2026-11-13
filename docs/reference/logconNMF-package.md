# Log-Concave Nonnegative Matrix Factorization with Automatic Rank Selection

Separates overlapping ordered signals into nonnegative log-concave
profiles and channel loadings, with automatic rank selection, optional
smooth backgrounds and flexible noise models. Estimates constant or
time-varying point spread functions and fits simultaneous models that
share loadings across samples on their original observation grids.

## Details

For a scans-by-channels signal, each component is the product of an
ordered nonnegative profile and a nonnegative vector of channel
loadings. The positive part of the profile has a concave logarithm,
which accommodates asymmetric peaks and tails. Chromatographic elution
profiles and mass spectra are a motivating example.

Use
[`estimatePSF`](https://tomwenseleers.github.io/logconNMF/reference/estimatePSF.md)
to build a sparse design from a standard or overlapping signals.
[`logconNMF`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)
estimates the number of components, their profiles and their loadings.
It fits one matrix or shares spectra across multiple samples with their
own profiles, areas and backgrounds.

## See also

Useful links:

- <https://tomwenseleers.github.io/logconNMF/>

- <https://github.com/tomwenseleers/logconNMF>

- Report bugs at <https://github.com/tomwenseleers/logconNMF/issues>

## Author

**Maintainer**: Tom Wenseleers <tom.wenseleers@kuleuven.be>

Authors:

- Tom Wenseleers <tom.wenseleers@kuleuven.be>

- Christophe Vanderaa <christophe.vanderaa@kuleuven.be>
