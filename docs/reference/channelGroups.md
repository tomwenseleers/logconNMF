# Learn and reuse a common channel aggregation

Cluster outcome channels using size-adjusted Poisson distances and Ward
linkage. The resulting map can be reused across samples or learned from
reference spectra supplied as rows with `normalize=TRUE`.

## Usage

``` r
channelGroups(
  y,
  k = 32L,
  channels = colnames(y),
  normalize = FALSE,
  ncores = NULL,
  reserve_unknown = TRUE
)

aggregateChannels(y, groups, channels = colnames(y))
```

## Arguments

- y:

  Numeric matrix with channels in columns.

- k:

  Maximum number of channel groups.

- channels:

  Unique channel identifiers, normally m/z values.

- normalize:

  Give every input row unit total, useful for spectral libraries.

- ncores:

  Native thread count, defaulting to physical cores.

- reserve_unknown:

  Reserve one group for channels absent from the training axis. Set to
  FALSE when fitting on a fixed, shared channel axis.

- groups:

  A map returned by `channelGroups()`.

## Value

`channelGroups()` returns a compact reusable map accepted as
`precluster` by
[`logconNMF()`](https://tomwenseleers.github.io/logconNMF/reference/logconNMF.md).
`aggregateChannels()` returns a matrix whose columns sum the mapped
channels. Unseen channels are collected in the map's reserved group.
Save maps with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

## Examples

``` r
y <- matrix(rpois(100*12, 3), 100, 12)
groups <- channelGroups(y, k=4L, ncores=1L)
z <- aggregateChannels(y, groups)
stopifnot(sum(y) == sum(z))
```
