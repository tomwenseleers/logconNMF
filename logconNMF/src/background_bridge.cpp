#include <Rcpp.h>
// [[Rcpp::depends(L0glm)]]
inline SEXP checked(SEXP result) {
  if (Rf_inherits(result,"try-error")) Rcpp::stop(Rcpp::as<std::string>(result));
  if (Rf_inherits(result,"interrupted-error")) throw Rcpp::internal::InterruptedException();
  if (Rcpp::internal::isLongjumpSentinel(result)) throw Rcpp::LongjumpException(result);
  return result;
}
// [[Rcpp::export]]
SEXP nuisance_poisson_cpp(SEXP B,SEXP Y,SEXP Signal,SEXP Prior,SEXP Start,
                          int maxit=60,double tol=1e-8,int nthreads=4) {
  typedef SEXP(*Fn)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
  static Fn fun=(Fn)R_GetCCallable("L0glm","_L0glm_nuisance_poisson_cpp");
  return checked(fun(B,Y,Signal,Prior,Start,Rcpp::wrap(maxit),Rcpp::wrap(tol),Rcpp::wrap(nthreads)));
}
// [[Rcpp::export]]
SEXP nuisance_poisson_checked_cpp(SEXP B,SEXP Y,SEXP Signal,SEXP Prior,SEXP Start,
                          int maxit=60,double tol=1e-8,int nthreads=4) {
  typedef SEXP(*Fn)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
  static Fn fun=(Fn)R_GetCCallable("L0glm","_L0glm_nuisance_poisson_checked_cpp");
  return checked(fun(B,Y,Signal,Prior,Start,Rcpp::wrap(maxit),Rcpp::wrap(tol),Rcpp::wrap(nthreads)));
}
// [[Rcpp::export]]
SEXP nuisance_gaussian_cpp(SEXP B,SEXP Y,SEXP Signal,SEXP Prior,SEXP Start,
                          SEXP constrained,int nthreads=4) {
  typedef SEXP(*Fn)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
  static Fn fun=(Fn)R_GetCCallable("L0glm","_L0glm_nuisance_gaussian_cpp");
  return checked(fun(B,Y,Signal,Prior,Start,constrained,Rcpp::wrap(nthreads)));
}
