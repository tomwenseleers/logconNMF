# Post-matching coordinate sweep on the original channels. The spectral and
# background updates see the same local rates and area normalization as before.
sinmf_shared_profile <- function(task, f, meta, old_threads, gid) {
  if (length(task$selected)) {
    cached <- sinmf_original_cache(task$directory, meta$n, meta$m)
    S <- task$P * (f$area * task$scale)
    if (is.null(f$profile_limits)) {
      padding <- ceiling(2 * f$small_width /
        (sqrt(log(100)/log(2)) * stats::median(diff(f$time))))
      f$profile_limits <- vapply(seq_along(gid), function(k) {
        at <- which(f$H[,k] > 0)
        if (!length(at)) return(c(1L,0L))
        as.integer(c(max(1,min(at)-padding), min(meta$n,max(at)+padding)))
      }, integer(2))
    }
    order <- task$selected[order(f$area[task$selected]*task$scale[task$selected], decreasing=TRUE)]
    for (k in order) {
      bounds <- f$profile_limits[,k]
      ions <- which(S[k,] > 0)
      if (!length(ions) || bounds[2L] < bounds[1L]) next
      at <- seq.int(bounds[1L], bounds[2L])
      h <- as.numeric(f$H[at,k])
      # Only overlapping components contribute to the local offset.
      local <- f$H[at,,drop=FALSE]
      neighbours <- which(Matrix::colSums(local)>0)
      neighbours <- neighbours[neighbours != k]
      off <- (f$background[at,,drop=FALSE] %*% f$alpha[,ions,drop=FALSE])*f$raw_exposure[at]
      if (length(neighbours))
        off <- off + as.matrix(local[,neighbours,drop=FALSE] %*% S[neighbours,ions,drop=FALSE])*f$exposure[at]
      slice <- function(x) if(is.null(dim(x))) sinmf_observation_block(x,at) else x[at,ions,drop=FALSE]
      off <- off + slice(f$offset)
      if (is.null(cached)) {
        y <- sinmf_read_matrix_file(file.path(task$directory,"observed.bin"),meta$n,meta$m,ions,at)
        columns <- seq_along(ions)
        first <- 0L
      } else {
        y <- cached
        columns <- ions
        first <- at[1L]-1L
      }
      rate <- slcn_fit_profile(y,columns,first,off,slice(f$weights),h,S[k,ions],
                 f$aligned_time[at],f$exposure[at],f$kind,sinmf_slice_noise(f$phi,ions),
                 sinmf_slice_noise(f$theta,ions),f$censor,task$extra$steps,f$ncores,TRUE)
      if (identical(rate,h)) next
      norm <- sum(rate*f$quadrature[at])
      if (!is.finite(norm) || norm<=0) next
      f$H[at,k] <- rate/norm
      f$area[k] <- f$area[k]*norm
      S[k,] <- task$P[k,]*f$area[k]*task$scale[k]
    }
  }
  f$gid <- gid
  f$ncores <- old_threads
  sinmf_checkpoint_sample(f,task$directory,"profile")
}
