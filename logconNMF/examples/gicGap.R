t <- seq_len(80)
H <- cbind(exp(-((t-25)/5)^2), exp(-((t-53)/4)^2))
S <- rbind(c(100,30,0), c(0,20,80))
Y <- round(2 + H %*% S)
X <- convMatrix(convKernel("gaussian",15L,list(sigma=2)),length(t))
fit <- logconNMF(Y,X,time=t,truth.profiles=H,ncores=1L,correlation="none")
fit$gic_gap$percent
plot(fit$gic_gap)
