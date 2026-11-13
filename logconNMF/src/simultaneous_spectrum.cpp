#include <RcppEigen.h>
#include <L0glm/scalar_families.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
// Pooled shared-spectrum score and Fisher information on the original scale.
// [[Rcpp::export]]
List slcn_pooled_terms(NumericMatrix y,NumericMatrix mean,NumericVector design,
 int family,NumericVector phi,NumericVector theta,NumericMatrix weights,int nthreads=1,double censor=.5) {
 const int n=y.nrow(),m=y.ncol();
 if(mean.nrow()!=n||mean.ncol()!=m||design.size()!=n||weights.nrow()!=n||!(weights.ncol()==1||weights.ncol()==m))
   stop("Incompatible pooled-spectrum dimensions.");
 if(!(phi.size()==1||phi.size()==m)||!(theta.size()==1||theta.size()==m)||
    nthreads<1||family<0||family>6)stop("Invalid pooled-spectrum parameters.");
 NumericVector gradient(m),information(m);
 #ifdef _OPENMP
 #pragma omp parallel for schedule(static) num_threads(nthreads)
 #endif
 for(int j=0;j<m;++j) {
   long double g=0,I=0;
   for(int i=0;i<n;++i)if(weights(i,weights.ncol()==1?0:j)>0&&design[i]>0) {
     auto z=L0glmScalar::derivatives(y(i,j),std::max(1e-15,mean(i,j)),family,phi[phi.size()==1?0:j],theta[theta.size()==1?0:j],censor);
     g+=weights(i,weights.ncol()==1?0:j)*design[i]*z.gradient;
     I+=weights(i,weights.ncol()==1?0:j)*design[i]*design[i]*z.information;
   }
   gradient[j]=(double)g;information[j]=(double)I;
 }
 return List::create(_["gradient"]=gradient,_["information"]=information);
}

// Pool all profiles in one pass per mass channel. Each scalar likelihood
// derivative is evaluated once per observed cell, independent of the rank.
// [[Rcpp::export]]
List slcn_pooled_spectra(NumericMatrix y,NumericMatrix mean,
 const Eigen::MappedSparseMatrix<double>& profiles,NumericVector exposure,
 NumericVector amplitude,int family,NumericVector phi,NumericVector theta,
 NumericVector weights,int nthreads=1,double censor=.5) {
 const int n=y.nrow(),m=y.ncol(),K=profiles.cols();
 if(mean.nrow()!=n||mean.ncol()!=m||profiles.rows()!=n||exposure.size()!=n||amplitude.size()!=K||
    !(weights.size()==1||weights.size()==n||weights.size()==n*m)||
    !(phi.size()==1||phi.size()==m)||!(theta.size()==1||theta.size()==m)||
    nthreads<1||family<0||family>6)stop("Invalid pooled profile dimensions.");
 NumericMatrix gradient(K,m),information(K,m);
 #ifdef _OPENMP
 #pragma omp parallel for schedule(static) num_threads(nthreads)
 #endif
 for(int j=0;j<m;++j) {
   std::vector<double> g(n),I(n),w(n);
   for(int i=0;i<n;++i) {
     w[i]=weights[weights.size()==1?0:i+(weights.size()==n?0:n*j)];
     if(w[i]<=0)continue;
     const auto z=L0glmScalar::derivatives(y(i,j),std::max(1e-15,mean(i,j)),family,
       phi[phi.size()==1?0:j],theta[theta.size()==1?0:j],censor);
     g[i]=z.gradient;I[i]=z.information;
   }
   for(int k=0;k<K;++k) {
     long double U=0,F=0;
     for(Eigen::MappedSparseMatrix<double>::InnerIterator it(profiles,k);it;++it) {
       const int i=it.row();const double x=(it.value()*exposure[i])*amplitude[k];
       if(w[i]>0&&x>0){U+=w[i]*x*g[i];F+=w[i]*x*x*I[i];}
     }
     gradient(k,j)=(double)U;information(k,j)=(double)F;
   }
 }
 return List::create(_["gradient"]=gradient,_["information"]=information);
}
