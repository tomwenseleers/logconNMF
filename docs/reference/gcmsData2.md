# GC/MS Chromatogram Section - Beer Headspace Volatiles

A 6 min segment (RT 19.00-25.00 min) of a GC/MS chromatogram recorded
from the SPME-headspace of beer. Data were acquired on a Thermo ISQ
single-quad coupled to a 30 m x 0.25 mm x 0.25 um RXi-5Sil MS column
using a DVB/Carboxen/PDMS (50/30 um) fibre. Many volatiles co-elute in
this interval, including by a prominent, triangular, "shark-fin"
octanoic-acid peak whose apex lies at RT 21.50 min. For details see
Schreurs et al. Nature Communications 2024.

## Usage

``` r
data(gcmsData2)
```

## Format

A numeric matrix with 1764 rows (time slices) and 518 columns (m/z
values).

Rows

: Individual full-scan acquisitions (time scans). Row names give the
  scan mid-point in minutes (19.00 \<= RT \<= 25.00).

Columns

: Unit-mass channels, m/z 33:550. Column names are the integer m/z
  values.

## Source

Schreurs, M. et al. (2024) "Predicting and improving complex beer flavor
through machine learning." \*Nature Communications\* \*\*15\*\*, 2368.
[doi:10.1038/s41467-024-46346-0](https://doi.org/10.1038/s41467-024-46346-0)

## Author

Tom Wenseleers <tom.wenseleers@kuleuven.be> \\ Laboratory of
Socioecology and Social Evolution, KU Leuven \\
<https://bio.kuleuven.be/eeb/tw/research> \\ ORCID:
[0000-0002-1434-861X](https://orcid.org/0000-0002-1434-861X)

## Examples

``` r
library(L0glm)
data(gcmsData2)
plotMat(gcmsData2,
      xlab = "Elution time (min)",
      ylab = "Mass (m/z)",
      main = "Section of beer headspace GC/MS profile")

matplot(x = as.numeric(rownames(gcmsData2)),
        y = gcmsData2^0.2, type = "l", lty = 1, col = rgb(0, 0, 0, 0.1),
        xlab = "Elution time (min)",
        ylab = "Intensity ^ 0.2",
        main = "Section of beer headspace GC/MS profile")

plot(x=as.numeric(names(attr(gcmsData2, "alkanes_TIC"))),
     y=attr(gcmsData2, "alkanes_TIC"),
     xlab="Retention time (min)", 
     ylab="Total ion intensity", 
     main="C7 to C22 alkane ladder standard", 
     type="l")
```
