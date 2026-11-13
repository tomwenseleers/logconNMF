# Run from the logconNMF repository root.
library(logconNMF)
RhpcBLASctl::blas_set_num_threads(1L)
results <- list()
for (name in c("gcmsData1","gcmsData2")) {
  data(list=name, package="logconNMF")
  Y <- get(name)
  time <- as.numeric(rownames(Y))
  X <- convMatrix(convKernel("gaussian",65L,
    list(sigma=if(name=="gcmsData1")8 else 6)),nrow(Y))
  B <- if(name=="gcmsData1") matrix(1,nrow(Y),1) else
    splines::bs(time,df=4L,intercept=TRUE)
  seconds <- system.time(fit <- splogconNMF(X,Y,background=B,
    family=quasipoisson("identity"),precluster=if(name=="gcmsData1")8L else 32L,
    correlation="auto", ncores=8L))[["elapsed"]]
  results[[name]] <- data.frame(dataset=name,n=nrow(Y),channels=ncol(Y),
    compounds=fit$rank,GIC=GIC(fit),seconds=seconds)
}
do.call(rbind,results)
