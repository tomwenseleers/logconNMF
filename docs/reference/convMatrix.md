# Make banded or rectangular convolution design matrices

`convMatrix()` builds a sparse Toeplitz convolution matrix \\A\\ such
that \$\$A x = \operatorname{conv}(x,g;\mathrm{shape}),\$\$ where \\g\\
is the supplied point-spread-function kernel.

- If `n_cols` is *omitted*, the function returns the classical *banded*
  design matrix: square when `shape = "same"` and tall when
  `shape = "full"`.

- If `n_cols` is supplied, only the first `n_cols` time-shifts of
  `kernel` are kept, yielding an \\n \times n\_{\mathrm{cols}}\\
  rectangular design. Such matrices are useful for point-spread-function
  (PSF) update steps in blind-deconvolution algorithms.

The matrix is produced in a single pass, so no intermediate (larger)
object is ever allocated, keeping memory usage to a minimum.

## Usage

``` r
convMatrix(kernel, n, shape = c("same", "full"), n_cols = NULL, origin = NULL)
```

## Arguments

- kernel:

  Numeric vector (\`length L\`) - FIR filter / PSF \*\*in forward-time
  order\*\*.

- n:

  Integer \\\>\\ 0 - length of the observed signal (the number of rows
  for a \*"same"\* design).

- shape:

  \`"same"\` (default, \\n\\ rows) or \`"full"\` (\\n+L-1\\ rows).

- n_cols:

  Optional. Number of columns to keep. \*If omitted\*,
  \\n\_{\mathrm{cols}}=n\\, and the result is the usual banded
  convolution matrix.

- origin:

  Integer (1 ... \\L\\) - element of \`kernel\` that is aligned with
  time 0. Defaults to the centre for \`"same"\` and to 1 for \`"full"\`.

## Value

A sparse matrix of class \*\*\`dgCMatrix\`\*\* with

\* \\n\_{\mathrm{rows}}=n\\ for shape \`"same"\`; or \*
\\n\_{\mathrm{rows}}=n+L-1\\ for shape \`"full"\`.

rows, and \`n_cols\` columns.

## See also

[`convKernel`](https://tomwenseleers.github.io/logconNMF/reference/convKernel.md)
to generate the 1-D PSF vectors used here.

## Author

Tom Wenseleers <tom.wenseleers@kuleuven.be> \\ Laboratory of
Socioecology and Social Evolution, KU Leuven \\
<https://bio.kuleuven.be/eeb/tw/research> \\ ORCID:
[0000-0002-1434-861X](https://orcid.org/0000-0002-1434-861X)

## Examples

``` r
## 1. synthetic spike-train convolved with a bi-Gaussian PSF ------------
if (requireNamespace("gsignal", quietly = TRUE)) {
  library(Matrix)
  library(gsignal)

  ### ---- synthetic spike train ----------------------------------------
  set.seed(1)
  n <- 100
  a <- numeric(n)
  a[sample(n, 20)] <- rpois(20, 10)      # random spikes

  ### ---- bi-Gaussian PSF via kernel() ---------------------------------
  L <- 31
  g <- convKernel("bigaussian", L,
                  params = list(sigma1 = 1, sigma2 = 4))  # default p = 0.5
  plot(g, type = "l", main = "Bi-Gaussian PSF", ylab = "density")

  ### ---------- 1.1 FULL convolution -----------------------------------
  y_ref_full <- as.vector(conv(a, g, shape = "full"))        # reference

  A_full <- convMatrix(g, n, shape = "full")                 # sparse full conv matrix
  all.equal(as.matrix(A_full), convmtx(g, n), tol = 1e-12)   # matches convmtx

  y_full <- as.vector(A_full %*% a)                          # A  %*%  a

  plot(1:(n + L - 1), y_ref_full, type = "l",
       main = "FULL convolution", xlab = "time index", ylab = "signal")
  lines(1:(n + L - 1), y_full, col = 2, lty = 3)             # identical

  centre_shift <- (L + 1) / 2 - 1                            # PSF peak offset
  lines((1:n) + centre_shift - 1, a, type = "h", col = 4)    # spike locations
  all.equal(y_full, y_ref_full, tol = 1e-12)                 # matches conv

  ### ---------- 1.2 SAME-size convolution ------------------------------
  y_ref_same <- as.vector(conv(a, g, shape = "same"))

  A_same <- convMatrix(g, n, shape = "same")                 # sparse conv matrix
  y_same <- as.vector(A_same %*% a)

  plot(1:n, y_ref_same, type = "l",
       main = "SAME-length convolution", xlab = "time index", ylab = "signal")
  lines(1:n, y_same, col = 2, lty = 3)                       # identical
  lines((1:n) - 1, a, type = "h", col = 4)                   # spike locations
  all.equal(y_same, y_ref_same, tol = 1e-12)                 # matches conv
}

## 2. rectangular design for PSF-update step in blinv deconv algo -------
set.seed(1)
n      <- 1000                                  # length of obs convolved signal
Lpsf   <- 101                                   # length of PSF to be estimated
spikes <- rpois(n, 2)                           # true or estimated spikes

G <- convMatrix(kernel = spikes,
                n      = n,
                n_cols = Lpsf,         # keep just the first Lpsf shifts
                shape  = "same",
                origin = (Lpsf + 1L) %/% 2L)
all.equal(            # identical to the first Lpsf columns of the banded design
  G,
  convMatrix(spikes, n, shape = "same",
             origin = (Lpsf + 1L) %/% 2L)[, 1:Lpsf, drop = FALSE], tol = 0)
```
