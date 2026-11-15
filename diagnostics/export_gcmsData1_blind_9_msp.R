#!/usr/bin/env Rscript

# Reproduce the exploratory nine-component blind-PSF fit discussed during the
# gcmsData1 chemical audit and export its resolved spectra for NIST MS Search.

pkgload::load_all("logconNMF", quiet = TRUE)

data(gcmsData1, package = "logconNMF")

# Current package data use the corrected nominal mass axis. Retain support for
# reproducing the export from a checkout predating that data correction.
if (identical(dim(gcmsData1), c(119L, 671L))) {
  stopifnot(
    identical(colnames(gcmsData1), as.character(30:700)),
    sum(gcmsData1[, 1L]) == 0
  )
  Y <- gcmsData1[, -1L, drop = FALSE]
  colnames(Y) <- as.character(as.integer(colnames(Y)) - 1L)
} else {
  stopifnot(
    identical(dim(gcmsData1), c(119L, 670L)),
    identical(colnames(gcmsData1), as.character(30:699))
  )
  Y <- gcmsData1
}
time <- as.numeric(rownames(Y))
mass <- as.integer(colnames(Y))

initial_kernel <- convKernel("gaussian", 41L, list(sigma = 8))
psf <- estimatePSF(
  Y,
  time = time,
  family = quasipoisson("identity"),
  standard = NULL,
  method = "blind",
  nwindows = 1L,
  precluster = 8L,
  ncores = 8L,
  control = list(
    L = 41L,
    kernel_est = initial_kernel,
    maxiter = 6L
  ),
  verbose = TRUE
)

fit <- logconNMF(
  Y,
  X = psf$X,
  background = matrix(1, nrow(Y), 1L),
  family = quasipoisson("identity"),
  precluster = 8L,
  correlation = "auto",
  ncores = 8L
)

if (fit$rank != 9L) {
  stop(
    "The exploratory fit did not reproduce rank nine (obtained rank ",
    fit$rank, ")."
  )
}

S <- as.matrix(coef(fit))
apex_scan <- if (!is.null(fit$apex_scan)) {
  as.integer(fit$apex_scan)
} else {
  max.col(t(as.matrix(fit$H)), ties.method = "first")
}
apex_rt <- time[apex_scan]
ion_threshold <- 0.001

output_file <- file.path(
  "diagnostics",
  "gcmsData1_blind_9_components_mass_corrected.msp"
)

entries <- unlist(lapply(seq_len(fit$rank), function(k) {
  spectrum <- S[k, ]
  keep <- which(spectrum >= ion_threshold * max(spectrum) & spectrum > 0)
  abundance <- 999 * spectrum[keep] / max(spectrum[keep])
  c(
    sprintf("Name: gcmsData1_blind_component_%02d_RT%.5f", k, apex_rt[k]),
    sprintf("Retention_time: %.7g", apex_rt[k]),
    sprintf("Pseudocount: %.15g", sum(spectrum)),
    sprintf("Ion_threshold: %.9g", ion_threshold),
    paste0(
      "Comments: Exploratory logconNMF blind-PSF component; identity unassigned; ",
      sprintf("apex retention time %.7g min; ", apex_rt[k]),
      sprintf("PSF median FWHM %.7g scans; ", psf$width_scans),
      "mass axis corrected to nominal m/z 30:699; ",
      "a factor is not by itself evidence for a distinct compound"
    ),
    sprintf("Num Peaks: %d", length(keep)),
    sprintf("%d %.8g;", mass[keep], abundance),
    ""
  )
}), use.names = FALSE)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
connection <- file(output_file, open = "wb")
on.exit(close(connection), add = TRUE)
writeLines(entries, connection, sep = "\r\n", useBytes = TRUE)
close(connection)
on.exit(NULL, add = FALSE)

summary <- data.frame(
  component = seq_len(fit$rank),
  apex_scan = apex_scan,
  apex_rt_min = apex_rt,
  exported_ions = vapply(seq_len(fit$rank), function(k) {
    sum(S[k, ] >= ion_threshold * max(S[k, ]) & S[k, ] > 0)
  }, integer(1L)),
  stringsAsFactors = FALSE
)

cat("Exported:", normalizePath(output_file, winslash = "/"), "\n")
cat("Rank:", fit$rank, "\n")
cat("PSF median FWHM (scans):", psf$width_scans, "\n")
print(summary, row.names = FALSE)

