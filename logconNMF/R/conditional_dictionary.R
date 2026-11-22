cc_width <- function(x) {
  if (!nrow(x$grid))
    return(0)
  median(x$warp[x$grid[, 2L]] - x$warp[x$grid[, 1L]])
}
cc_unit <- function(p) p/pmax(sqrt(rowSums(p * p)), 1e-80)
cc_area_shape <- function(t, h) {
  if (!any(h > 0))
    return(h)
  z <- cc_log_shape(t, h)
  z * (sum(h)/max(sum(z), 1e-80))
}
cc_means <- function(x, p) {
  # A chromatographic profile occupies a small interval of a long scan grid.
  x$mu <- pmax(x$bg + as.matrix(as(x$h, "dgCMatrix") %*% p), 1e-08)
}
cc_lift <- function(x, p, group, price, threads) {
  count <- rowSums(p)
  x$h <- as.matrix(x$init %*% Matrix::Diagonal(x = count))
  p <- p/pmax(count, 1e-100)
  design <- cbind(as(x$h, "dgCMatrix"), as(x$z, "dgCMatrix"))
  w <- 1/sqrt(x$phi * pmax(rowSums(x$y), 1))
  X <- Matrix::Diagonal(x = w) %*% design
  scale <- sqrt(Matrix::colSums(X^2))
  scale[scale == 0] <- 1
  X <- X %*% Matrix::Diagonal(x = 1/scale)
  ion <- colSums(x$y)
  denom <- numeric(max(group))
  sums <- rowsum(ion, group, reorder = FALSE)
  denom[as.integer(rownames(sums))] <- sums[, 1]
  frac <- ion/pmax(denom[group], 1e-100)
  initial <- sweep(rbind(p, x$alpha)[, group, drop = FALSE], 2L, frac, "*")
  coef <- cc_gram_nnls(as(Matrix::crossprod(X), "generalMatrix"), as.matrix(Matrix::crossprod(X, x$y * w)), initial *
    scale)/scale
  K <- nrow(p)
  idx <- seq_len(K)
  p <- coef[idx, , drop = FALSE]
  size <- rowSums(p)
  p <- p/pmax(size, 1e-100)
  x$h <- sweep(x$h, 2, size, "*")
  x$alpha <- coef[K + seq_len(ncol(x$z)), , drop = FALSE]
  x$bg <- pmax(x$z %*% x$alpha + 1e-08, 1e-08)
  cc_means(x, p)
  for (i in 1:2) cc_spectrum_sweep(list(x), p, order(colSums(x$h), decreasing = TRUE), price, FALSE,
    threads)
  area <- colSums(x$h)
  live <- which(area > 0)
  x$initial_p <- p[live, , drop = FALSE] * area[live]
  x$grid <- x$grid[live, , drop = FALSE]
  x$init <- as(x$h[, live, drop = FALSE], "dgCMatrix") %*% Matrix::Diagonal(x = 1/area[live])
  x$h <- x$mu <- NULL
}
cc_sample <- function(part, group, price, threads) {
  x <- new.env(parent = emptyenv())
  x$init <- as(part$H, "dgCMatrix")
  x$initial_p <- as.matrix(part$S)
  x$y <- if (is.character(part$observed))
    sinmf_read_matrix_file(part$observed, nrow(x$init), length(group)) else as.matrix(part$observed)
  x$z <- unclass(as.matrix(part$background))
  if (!is.double(x$y)) storage.mode(x$y) <- "double"
  if (!is.double(x$z)) storage.mode(x$z) <- "double"
  x$alpha <- as.matrix(part$alpha) + 0
  x$t <- as.double(part$time)
  x$warp <- as.double(part$aligned)
  x$grid <- if (is.null(part$windows)) conditional_windows(x$init) else part$windows
  x$phi <- part$spec$dispersion
  n <- length(x$t)
  if (!identical(dim(x$y), c(n, length(group))) || length(x$warp) != n ||
      any(!is.finite(x$t)) || any(!is.finite(x$warp)) || any(diff(x$t) <= 0) ||
      any(diff(x$warp) <= 0) || nrow(x$init) != n || nrow(x$z) != n ||
      ncol(x$z) != nrow(x$alpha) || nrow(x$initial_p) != ncol(x$init) ||
      ncol(x$initial_p) != ncol(x$alpha) || !is.finite(x$phi) || x$phi <= 0)
    stop("Invalid count dictionary or registration.")
  x$width <- max(cc_width(x), if (length(x$t) > 1) 3 * median(diff(x$warp)) else 0)
  if (isTRUE(part$grouped_dictionary) || ncol(x$initial_p) != ncol(x$y))
    cc_lift(x, x$initial_p, group, price, threads) else x$bg <- pmax(x$z %*% x$alpha + 1e-08, 1e-08)
  x
}
cc_atoms <- function(data) {
  a <- lapply(seq_along(data), function(i) {
    x <- data[[i]]
    g <- x$grid
    x$norms <- rowSums(x$initial_p)
    heights <- vapply(seq_len(ncol(x$init)), function(j) {
      pos <- seq.int(x$init@p[j] + 1L, x$init@p[j + 1L])
      max(x$init@x[pos])
    }, 0)
    data.frame(si = rep.int(i, nrow(g)), id = seq_len(nrow(g)), lo = x$warp[g[, 1]], hi = x$warp[g[,
      2]], peak = x$warp[g[, 3]], height = heights * x$norms)
  })
  list(a = do.call(rbind, a), p = do.call(rbind, lapply(data, function(x) x$initial_p/pmax(x$norms,
    1e-80))))
}
cc_dictionary <- function(data, price, threads) {
  atoms <- cc_atoms(data)
  a <- atoms$a
  p <- atoms$p
  K <- nrow(p)
  m <- ncol(p)
  bands <- as.matrix(a[, c("lo", "hi")])
  dimnames(bands) <- NULL
  if (!K) {
    for (x in data) x$h <- matrix(0, length(x$t), 0)
    return(list(p = p, bands = bands))
  }
  if (length(data) == 1L) {
    x <- data[[1L]]
    x$h <- as.matrix(x$init %*% Matrix::Diagonal(x = x$norms))
    for (j in seq_len(K)) {
      # Preserve an exact feasible starting shape. Projection of feasible curves is the identity.
      if (!cc_feasible(x$t, x$h[, j]))
        x$h[, j] <- cc_area_shape(x$t, x$h[, j])
    }
    if (m == 1L)
      return(list(p = p, bands = bands))
    unit <- cc_unit(p)
    todo <- rep(TRUE, K)
    gr <- list()
    for (j in order(-a$height)) if (todo[j]) {
      ok <- todo & drop(unit %*% unit[j, ]) >= 0.97
      lo <- a$lo[j]
      hi <- a$hi[j]
      repeat {
        ids <- which(ok & a$hi >= lo & a$lo <= hi)
        left <- min(a$lo[ids])
        right <- max(a$hi[ids])
        if (left == lo && right == hi)
          break
        lo <- left
        hi <- right
      }
      todo[ids] <- FALSE
      gr[[length(gr) + 1L]] <- ids
    }
    if (length(gr) < K) {
      new <- matrix(0, length(gr), m)
      H <- matrix(0, length(x$t), length(gr))
      b <- matrix(0, length(gr), 2L)
      for (i in seq_along(gr)) {
        g <- gr[[i]]
        w <- sqrt(a$height[g])
        new[i, ] <- colSums(p[g, , drop = FALSE] * w)/sum(w)
        H[, i] <- cc_area_shape(x$t, as.numeric(x$init[, g, drop = FALSE] %*% x$norms[g]))
        b[i, ] <- c(min(a$lo[g]), max(a$hi[g]))
      }
      cost <- function(h, q) {
        mu <- x$bg + as.matrix(as(h, "dgCMatrix") %*% q)
        cc_deviance(x$y, mu, x$phi) + price * (sum(cc_faces(x$t, h)) + sum(q > 0) - nrow(q))
      }
      if (cost(H, new) < cost(x$h, p)) {
        x$h <- H
        return(list(p = new, bands = b))
      }
    }
    return(list(p = p, bands = bands))
  }
  raw <- p
  purity <- numeric(K)
  for (i in seq_len(K)) {
    x <- data[[a$si[i]]]
    g <- x$grid[a$id[i], ]
    core <- max(1, round((g[2] - g[1])/6))
    ix <- seq.int(max(1, g[3] - core), min(length(x$t), g[3] + core))
    y <- x$y[ix, , drop = FALSE]
    o <- x$bg[ix, , drop = FALSE]
    r <- pmax(y - o, 0)
    h <- rowSums(r)
    z <- cc_local_spectrum(y/x$phi, o/x$phi, h/x$phi, p[i, ], price)
    if (sum(z) > 0)
      raw[i, ] <- z/sum(z)
    err <- sum((r - tcrossprod(h, raw[i, ]))^2) + x$phi * sum(y)
    purity[i] <- sum(r * r)/max(err, 1)
  }
  src <- integer(K)
  u <- cc_unit(raw)
  b <- matrix(0, K, 2)
  n <- 0L
  for (i in order(-purity)) {
    ids <- seq_len(n)
    c <- ids[b[ids, 1] <= a$hi[i] & b[ids, 2] >= a$lo[i]]
    if (length(c)) {
      cos <- drop(u[src[c], , drop = FALSE] %*% u[i, ])
      c <- c[utils::head(order(-cos), 8L)]
      design <- t(u[src[c], , drop = FALSE])
      w <- cc_cone(design, u[i, ])
      if (sum((u[i, ] - design %*% w)^2) < 1 - 0.98^2) {
        j <- c[which.max(w)]
        if (sum(u[src[j], ] * u[i, ]) > 0.98)
          b[j, ] <- c(min(b[j, 1], a$lo[i]), max(b[j, 2], a$hi[i]))
        next
      }
    }
    n <- n + 1L
    src[n] <- i
    b[n, ] <- c(a$lo[i], a$hi[i])
  }
  ids <- seq_len(n)
  b <- b[ids, , drop = FALSE]
  q <- raw[src[ids], , drop = FALSE]
  u <- cc_unit(q)
  for (si in seq_along(data)) {
    x <- data[[si]]
    co <- matrix(0, nrow(x$initial_p), n)
    for (i in which(a$si == si)) {
      c <- which(b[, 1] <= a$hi[i] & b[, 2] >= a$lo[i])
      if (!length(c))
        next
      c <- c[utils::head(order(-drop(u[c, , drop = FALSE] %*% p[i, ])), 8L)]
      w <- cc_cone(t(u[c, , drop = FALSE]), p[i, ]/max(sqrt(sum(p[i, ]^2)), 1e-60)) * rowSums(u[c,
        , drop = FALSE])
      if (sum(w) > 0)
        co[a$id[i], c] <- x$norms[a$id[i]] * w/sum(w)
    }
    h <- as.matrix(x$init %*% co)
    mask <- matrix(FALSE, nrow(h), n)
    width <- cc_width(x)
    for (j in which(colSums(h) > 0)) {
      ix <- which(h[, j] >= 0.01 * max(h[, j]))
      left <- x$warp[ix[1]] - width
      right <- x$warp[utils::tail(ix, 1)] + width
      mask[, j] <- x$warp >= left & x$warp < right
    }
    h[!mask] <- 0
    h <- cc_scan_unmix(x$y, x$bg, q, h, mask, x$phi)
    for (j in which(colSums(h) > 0)) h[, j] <- cc_area_shape(x$t, h[, j])
    x$h <- h
  }
  list(p = q, bands = b)
}
cc_feasible <- function(t, h) {
  pos <- which(h > 0)
  if (length(pos) < 2L)
    return(TRUE)
  if (any(diff(pos) != 1L))
    return(FALSE)
  if (length(pos) < 3L)
    return(TRUE)
  z <- log(h[pos])
  t <- t[pos]
  j <- seq.int(2L, length(z) - 1L)
  all(z[j - 1L] + (z[j + 1L] - z[j - 1L]) * (t[j] - t[j - 1L])/(t[j + 1L] - t[j - 1L]) - z[j] <= 1e-09)
}
