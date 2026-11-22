# Run from the logconNMF repository root. Timings cover factorization with a
# fixed, boundary-aware PSF; CDF input/export is measured in GCMSprime examples.
library(logconNMF)
RhpcBLASctl::blas_set_num_threads(1L)
results <- list()
for (name in c("gcmsData1", "gcmsData2")) {
  data(list=name, package="logconNMF")
  Y <- get(name);time <- as.numeric(rownames(Y))
  first <- name=="gcmsData1"
  psf <- estimatePSF(Y, family=quasipoisson("identity"), time=time,
    standard=NULL, nwindows=1L, precluster=1L, ncores=8L,
    control=list(kernel_est=convKernel("gaussian",65L,
      list(sigma=if(first)9 else 6)),maxiter=0L))
  B <- if(first) matrix(1,nrow(Y),1L) else
    splines::bs(time,df=6L,intercept=TRUE)
  seconds <- system.time(fit <- logconNMF(Y,X=psf$X,background=B,time=time,
    family=quasipoisson("identity"),precluster=if(first)7L else 32L,
    correlation="auto",ncores=8L,
    control=list(initial_relaxed=!first,fusion_passes=0L)),gcFirst=FALSE)[["elapsed"]]
  results[[name]] <- data.frame(dataset=name,n=nrow(Y),channels=ncol(Y),
    compounds=fit$rank,GIC=GIC(fit),seconds=seconds)
}
result <- do.call(rbind,results)
dir.create("benchmarks/results",showWarnings=FALSE)
write.csv(result,"benchmarks/results/gcms_windows.csv",row.names=FALSE)
print(result)
