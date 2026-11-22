#include <Rcpp.h>
#include "count_projection.h"
#include <L0glm/scalar_families.h>
#include <algorithm>
#include <numeric>
#include <vector>

// Weighted projection of log rates onto the concave, piecewise-linear cone.
// With fixed knots the normal equations are tridiagonal.
namespace {
using namespace logcon_geometry;
bool feasible(const Vec& h, const Vec& t) {
    bool started=false, ended=false;
    for (size_t i=0;i<h.size();++i) {
        if (h[i]==0) { if (started) ended=true; continue; }
        if (ended) return false;
        started=true;
        if (i<2 || h[i-2]==0 || h[i-1]==0) continue;
        double b=(t[i-1]-t[i-2])/(t[i]-t[i-2]);
        if (std::log(h[i-1])+1e-10 < (1-b)*std::log(h[i-2])+b*std::log(h[i])) return false;
    }
    return true;
}

// Pointers to one ion's contiguous scan series. Scalar and scan-specific
// weights/offsets are broadcast by stride, without expanding an n x m matrix.
struct Ion {
    const double *y, *offset, *weight;
    int os, ws;
    double p, phi, theta;
};
}

// [[Rcpp::export]]
Rcpp::NumericVector slcn_fit_profile(Rcpp::NumericMatrix y, Rcpp::IntegerVector channels,
    int first, Rcpp::NumericVector offset, Rcpp::NumericVector weights,
    Rcpp::NumericVector start, Rcpp::NumericVector spectrum, Rcpp::NumericVector time,
    Rcpp::NumericVector exposure, int family, Rcpp::NumericVector phi,
    Rcpp::NumericVector theta, double censor, int maxit, int nthreads, bool require_improvement=false) {
    const int n=start.size(), m=spectrum.size();
    const auto ny=y.nrow();
    auto dimensions=[&](R_xlen_t len){ return len==1 || len==n || len==static_cast<R_xlen_t>(n)*m; };
    if (!n || !m || channels.size()!=m || first<0 || first>ny-n || time.size()!=n ||
        exposure.size()!=n || !dimensions(offset.size()) || !dimensions(weights.size()) ||
        !(phi.size()==1 || phi.size()==m) || !(theta.size()==1 || theta.size()==m) ||
        family<0 || family>6 || maxit<0 || nthreads<1) Rcpp::stop("Invalid profile dimensions or controls.");
    const double mass=std::accumulate(spectrum.begin(),spectrum.end(),0.);
    if (!std::isfinite(mass) || mass<=0) Rcpp::stop("Spectrum must have positive finite mass.");
    std::vector<Ion> ions;
    for (int j=0;j<m;++j) {
        if (channels[j]<1 || channels[j]>y.ncol()) Rcpp::stop("Invalid channel index.");
        if (spectrum[j]==0) continue;
        ions.push_back({y.begin()+static_cast<R_xlen_t>(channels[j]-1)*ny+first,
            offset.begin()+(offset.size()>n ? static_cast<R_xlen_t>(j)*n : 0),
            weights.begin()+(weights.size()>n ? static_cast<R_xlen_t>(j)*n : 0),
            offset.size()!=1,weights.size()!=1,spectrum[j]/mass,
            phi[phi.size()==1?0:j],theta[theta.size()==1?0:j]});
    }
    const bool poisson=(family==1 || family==2);
    const double phi_scale=poisson ? ions.front().phi : 1.;
    const double *e=exposure.begin();
    // R's special-function paths for censored Gamma and NB1 stay on one thread.
    const int threads=(family<=3 && n*double(ions.size())>=8192) ? std::min(n,nthreads) : 1;
    Vec t(n), old(n), free_rate(n), scratch(n);
    for (int i=0;i<n;++i) {
        t[i]=n==1?0:(time[i]-time[0])/(time[n-1]-time[0]);
        old[i]=start[i]*mass;
    }
    // Independent scan-wise starts. For Poisson, bracket the scalar likelihood
    // root. The other families use positive weighted least squares as a pilot.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
#endif
    for (int i=0;i<n;++i) {
        double numerator=0, denominator=0, observed_signal=0, sum_y=0, sum_p=0, upper=0;
        bool unit=true;
        for (const auto& ion:ions) {
            double w=ion.weight[i*ion.ws]/ion.phi, obs=ion.y[i], bg=ion.offset[i*ion.os];
            numerator+=w*ion.p*(obs-bg); denominator+=w*ion.p*ion.p;
            if (family>=3) observed_signal+=w*ion.p*obs;
            sum_y+=obs; sum_p+=w*ion.p;
            if (w>0) upper=std::max(upper,std::max(0.,obs-bg)/ion.p);
            unit=unit && w==ions.front().weight[i*ions.front().ws]/ions.front().phi;
        }
        if (!poisson || sum_p<=0) {
            // For variable-shape NB1 and Gamma the mean MLE can exceed the
            // observation even when it equals the background. A positive seed
            // lets their likelihoods compete with the exact zero candidate.
            double seed=numerator>0 ? numerator : .05*observed_signal;
            free_rate[i]=denominator>0?std::max(0.,seed/denominator)/e[i]:old[i];
            continue;
        }
        double lo=0, hi=unit?sum_y:upper, x=.25*hi;
        for (int k=0;k<26 && hi>0;++k) {
            double g=sum_p, curvature=0;
            for (const auto& ion:ions) {
                double w=ion.weight[i*ion.ws]/ion.phi;
                if (ion.y[i]==0 || w==0) continue;
                double r=ion.p/std::max(1e-150,ion.offset[i*ion.os]+x*ion.p);
                g-=w*ion.y[i]*r;
                curvature+=w*ion.y[i]*r*r;
            }
            if (std::abs(g)<1e-8*sum_p) break;
            if (g>=0) hi=x; else lo=x;
            double next=x-g/std::max(curvature,1e-100);
            x=(lo<next && next<hi)?next:.5*(lo+hi);
        }
        free_rate[i]=x/e[i];
    }
    // A common finite reference handles a zero warm start at positive counts.
    auto loss=[&](const Vec& h) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
#endif
        for (int i=0;i<n;++i) {
            double sum=0;
            for (const auto& ion:ions) {
                double w=ion.weight[i*ion.ws];
                if (w==0) continue;
                double yij=ion.y[i], bg=ion.offset[i*ion.os];
                double mu=bg+e[i]*h[i]*ion.p, ref=bg+e[i]*old[i]*ion.p;
                if (ref==0 && yij>0 && family!=0) ref=std::max(1.,yij);
                if (poisson) {
                    if (mu==0 && yij>0) {sum=R_PosInf;break;}
                    sum+=w*(phi_scale/ion.phi)*((mu-ref)-(yij>0?yij*(std::log(mu)-std::log(ref)):0));
                } else {
                    sum+=.5*w*L0glmScalar::change(yij,ref,mu,family,ion.phi,ion.theta,censor);
                }
            }
            scratch[i]=sum;
        }
        return std::accumulate(scratch.begin(),scratch.end(),0.);
    };
    Vec best(n,0);
    double best_loss=loss(best);
    if (feasible(old,t)) {
        double value=loss(old);
        if (value<=best_loss) {best=old;best_loss=value;}
    }
    for (const auto* seed : {&old,&free_rate}) {
        double peak=largest(*seed);
        if (!(peak>0)) continue;
        Vec u(n), w(n), rate(n), gradient(n), info(n), target(n);
        for (int i=0;i<n;++i) {
            u[i]=std::log(std::max((*seed)[i],peak*1e-9));
            w[i]=std::max((*seed)[i]/peak,1e-6);
        }
        u=project(t,u,w);
        for (int it=0;it<=maxit;++it) {
            for (int i=0;i<n;++i) rate[i]=std::exp(u[i]);
            double value=loss(rate);
            if (std::isfinite(value) && value<best_loss) {best_loss=value;best=rate;}
            if (it==maxit) break;
            Rcpp::checkUserInterrupt();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
#endif
            for (int i=0;i<n;++i) {
                double g=0, fisher=0;
                for (const auto& ion:ions) {
                    double wobs=ion.weight[i*ion.ws];
                    if (wobs==0) continue;
                    double signal=e[i]*rate[i]*ion.p;
                    double mu=std::max(1e-150,ion.offset[i*ion.os]+signal);
                    if (poisson) {
                        double ratio=signal/mu, scale=wobs*(phi_scale/ion.phi);
                        g+=scale*(signal-ion.y[i]*ratio); fisher+=scale*signal*ratio;
                    } else {
                        auto a=L0glmScalar::derivatives(ion.y[i],mu,family,ion.phi,ion.theta,censor);
                        g+=wobs*signal*a.gradient; fisher+=wobs*signal*signal*a.information;
                    }
                }
                gradient[i]=g;info[i]=fisher;
            }
            double floor=1e-9*std::max(1.,largest(info));
            for (int i=0;i<n;++i) {
                w[i]=std::max({floor,info[i],std::abs(gradient[i])/3.});
                target[i]=u[i]-gradient[i]/w[i];
            }
            Vec dir=project(t,target,w);
            for (int i=0;i<n;++i) dir[i]-=u[i];
            double slope=std::inner_product(gradient.begin(),gradient.end(),dir.begin(),0.);
            if (!(slope < -1e-9)) break;
            bool accepted=false;
            for (int half=0;half<16;++half) {
                double step=std::ldexp(1.,-half);
                for (int i=0;i<n;++i) rate[i]=std::exp(u[i]+step*dir[i]);
                if (loss(rate)<=value+1e-4*step*slope) {
                    for (int i=0;i<n;++i) u[i]+=step*dir[i];
                    accepted=true;break;
                }
            }
            if (!accepted) break;
        }
    }
    if (require_improvement && best_loss>loss(old)) return start;
    Rcpp::NumericVector answer(n);
    for (int i=0;i<n;++i) answer[i]=best[i]/mass;
    return answer;
}
