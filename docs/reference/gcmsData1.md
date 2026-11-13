# GC/MS Chromatogram Section - Cuticular Hydrocarbon Profile of a Cataglyphis iberica Ant Queen

A small 0.4 min section of a GC/MS chromatogram (RT 26.80-27.20 min)
obtained from the cuticular hydrocarbon extract of a \*Cataglyphis
iberica\* queen. Six compounds with near-Gaussian peak shapes co-elute
in this segment. Data were collected on a Shimadzu QP-2010 Ultra
equipped with a 30 m x 0.25 mm x 0.25 um DB-5ms column after liquid
injection of a hexane extract. For details see Van Oystaeyen et al.
Science 2014.

## Usage

``` r
data(gcmsData1)
```

## Format

A numeric matrix with 119 rows (time slices) and 670 columns (m/z
values).

Rows

: Individual full-scan acquisitions (time scans). Row names give the
  scan mid-point in minutes (26.80 \<= RT \<= 27.20).

Columns

: Unit-mass channels, m/z 30:699. Column names are the integer m/z
  values.

## Source

Van Oystaeyen, A. et al. (2014) "Conserved class of queen pheromones
stops social insect workers from reproducing." \*Science\* \*\*343\*\*,
287-290.
[doi:10.1126/science.1244899](https://doi.org/10.1126/science.1244899)

## Details

Co-eluting compounds present with retention times and diagnostic ions:

|  |  |  |
|----|----|----|
| **RT (min)** | **Compound** | **Diagnostic ions** |
| 26.84 | 15-Methylnonacosane | 224, 225 |
| 26.91 | 13-Methylnonacosane | 196, 252 |
| 26.95 | 11-Methylnonacosane | 168, 280 |
| 26.95 | Cyclodecasiloxane, eicosamethyl | 73, 147, 221, 281, 355, 429, 503 |
| 27.01 | 9-Methylnonacosane | 140, 308 |
| 27.11 | 7-Methylnonacosane | 112, 336 |

## Author

Tom Wenseleers <tom.wenseleers@kuleuven.be> \\ Laboratory of
Socioecology and Social Evolution, KU Leuven \\
<https://bio.kuleuven.be/eeb/tw/research> \\ ORCID:
[0000-0002-1434-861X](https://orcid.org/0000-0002-1434-861X)

## Examples

``` r
library(L0glm)
data(gcmsData1)
plotMat(gcmsData1,
      xlab = "Elution time (min)",
      ylab = "Mass (m/z)",
      main = "Section of Cataglyphis iberica\nant queen cuticular GC/MS profile")
matplot(x = as.numeric(rownames(gcmsData1)),
        y = gcmsData1^0.2, type = "l", lty = 1, col = rgb(0, 0, 0, 0.1),
        xlab = "Elution time (min)",
        ylab = "Intensity ^ 0.2",
        main = "Section of Cataglyphis iberica\nant queen cuticular GC/MS profile")
plot(x=as.numeric(names(attr(gcmsData1, "alkanes_TIC"))),
     y=attr(gcmsData1, "alkanes_TIC"),
     xlab="Retention time (min)", 
     ylab="Total ion intensity", 
     main="C9 to C40 alkane ladder standard", 
     type="l")
```
