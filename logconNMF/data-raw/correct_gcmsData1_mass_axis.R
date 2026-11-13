# Correct the historical one-unit mass-axis offset in gcmsData1.
#
# The stored matrix had an all-zero leading column labelled m/z 30, followed by
# intensities for nominal m/z 30:699 labelled as 31:700.  This script performs
# the minimal correction: discard the empty column, decrement the remaining
# labels, and preserve the alkane-standard TIC attribute verbatim.

data_file <- file.path("logconNMF", "data", "gcmsData1.rda")
load(data_file)

stopifnot(
  is.matrix(gcmsData1),
  identical(dim(gcmsData1), c(119L, 671L)),
  identical(colnames(gcmsData1), as.character(30:700)),
  all(gcmsData1[, 1L] == 0)
)

alkanes_TIC <- attr(gcmsData1, "alkanes_TIC", exact = TRUE)
gcmsData1 <- gcmsData1[, -1L, drop = FALSE]
colnames(gcmsData1) <- as.character(30:699)
attr(gcmsData1, "alkanes_TIC") <- alkanes_TIC

stopifnot(
  identical(dim(gcmsData1), c(119L, 670L)),
  identical(colnames(gcmsData1), as.character(30:699)),
  identical(attr(gcmsData1, "alkanes_TIC", exact = TRUE), alkanes_TIC)
)

save(gcmsData1, file = data_file, compress = "xz", version = 3L)
