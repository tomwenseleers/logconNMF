# Empirical GIC gap from matched reference profiles

The second fit starts from unit-area profiles on the observed scan grid.
A least concave majorant of log height gives a feasible starting shape
when irregular acquisition or bin integration alters log-concavity. Both
fits estimate spectra, amplitudes, shapes and background from the same
data. The block traces include the public GIC at each checkpoint. With a
different selection penalty, its decrease need not imply a GIC decrease.

## Usage

``` r
gicGap(object)

# S3 method for class 'logconNMF_gap'
plot(x, type = c("both", "GIC", "gap"), ...)

# S3 method for class 'logconNMF_gap'
print(x, ...)
```

## Arguments

- object:

  A logconNMF or GCMSprime fit computed with `truth.profiles`.

- x:

  The returned diagnostic.

- type:

  Show both panels, the complete GIC trace, or the final half of each
  trajectory on a percentage scale relative to the reference endpoint.

- ...:

  Other arguments passed to the base plot.

## Value

A plottable object with the standard and reference GIC, their signed
`difference`, `percent = 100 * difference / abs(reference_gic)`, both
block traces and elapsed seconds for reference initialization and
fitting. A zero reference GIC gives an undefined percentage. A negative
difference means the standard fit reached a lower GIC. The reference is
a feasible fitted solution. A global optimality bound would require a
separate proof. An independently fitted GCMSprime directory returns a
named list of these diagnostics, one per sample.

## See also

[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md)

## Examples

``` r
t <- seq_len(80)
H <- cbind(exp(-((t-25)/5)^2), exp(-((t-53)/4)^2))
S <- rbind(c(100,30,0), c(0,20,80))
Y <- round(2 + H %*% S)
X <- convMatrix(convKernel("gaussian",15L,list(sigma=2)),length(t))
fit <- logconNMF(Y,X,time=t,truth.profiles=H,ncores=1L,correlation="none")
fit$gic_gap$percent
#> [1] -40.38993
plot(fit$gic_gap)
```
