# Source-mode testthat runs always load setup files; define the shared numerical
# fixture here so filtered and full test runs use the same environment.
polish_example <- function(kind = 2L, distinct = FALSE, precision = 1) {
  n <- 120L;m <- 4L
  profile <- function(center, sigma, left, right) {
    h <- numeric(n);at <- seq.int(left,right)
    h[at] <- exp(-(at-center)^2/(2*sigma^2));h/sum(h)
  }
  H <- Matrix::Matrix(cbind(profile(65,18,15,110),
    profile(if(distinct)55 else 65,if(distinct)4 else 18,if(distinct)40 else 15,
      if(distinct)70 else 110)),sparse=TRUE)
  S <- rbind(c(200000,100000,10000,1000),
    if(distinct)c(20000,10000,1500,2000)else c(10000,5000,500,50))
  B <- matrix(1,n,1L);A <- matrix(2,1L,m);O <- matrix(.2,n,m)
  mu <- as.matrix(H%*%S)+B%*%A+O
  W <- matrix(precision,n,m);W[1:3,] <- 0;W[,2] <- W[,2]/2
  phi <- c(1,1,3,1,.5,3,3)[kind+1L]
  f <- list(H=H,S=S,background=B,alpha=A,offset=O,mu=mu,weights=W,
    spec=list(kind=kind,dispersion=phi,theta=20),lambda=1000,
    ncores=1L,control=list(censor=.5,psf_width=2L),members=list(1L,2L))
  Y <- if(kind%in%c(1L,3L,5L))round(mu)else mu
  list(f=slcn_sync(f,Y),Y=Y)
}

