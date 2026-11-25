#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
struct Sample {
    int scans;
    double dispersion;
    const double *observed, *mean, *profile;
};

// Half the change in Poisson deviance. Use log1p near the incumbent.
inline double loss_change(double y, double mu, double change) {
    const double next = std::max(1e-8, mu + change);
    const double delta = next - mu;
    const double log_ratio = std::abs(delta) < .5 * mu ?
        std::log1p(delta / mu) : std::log(next) - std::log(mu);
    return delta - y * log_ratio;
}
}

// Original-ion costs for deleting one occurrence, or transferring its profile
// to a neighbour or a convex mixture of two neighbours. A transferred profile
// is a screening proposal. The subsequent refit enforces the shape constraint.
// [[Rcpp::export]]
Rcpp::NumericMatrix cc_screen_transfers(Rcpp::List data, Rcpp::NumericMatrix spectra,
        Rcpp::IntegerMatrix first, Rcpp::IntegerMatrix second, Rcpp::NumericMatrix fraction,
        int threads) {
    const int rank = spectra.nrow(), ions = spectra.ncol(), samples = data.size();
    if (threads < 1 || first.nrow() != rank || first.ncol() != samples ||
        second.nrow() != rank || second.ncol() != samples ||
        fraction.nrow() != rank || fraction.ncol() != samples)
        Rcpp::stop("Invalid transfer dimensions or thread count.");
    for (int i : first) if (i < 0 || i > rank) Rcpp::stop("Invalid first recipient.");
    for (int i : second) if (i < 0 || i > rank) Rcpp::stop("Invalid second recipient.");
    for (double w : fraction) if (!R_FINITE(w) || w < 0 || w > 1)
        Rcpp::stop("Invalid transfer fraction.");
    std::vector<Sample> input;
    for (SEXP item : data) {
        Rcpp::Environment x(item);
        Rcpp::NumericMatrix h = x["h"], y = x["y"], mu = x["mu"];
        const double phi = Rcpp::as<double>(x["phi"]);
        if (h.ncol() != rank || y.ncol() != ions || mu.ncol() != ions ||
            h.nrow() != y.nrow() || y.nrow() != mu.nrow() || !R_FINITE(phi) || phi <= 0)
            Rcpp::stop("Incompatible count fit in support search.");
        input.push_back({h.nrow(), phi, REAL(y), REAL(mu), REAL(h)});
    }
    Rcpp::NumericMatrix result(rank, samples);
    const double *p = REAL(spectra), *weights = REAL(fraction);
    const int *one = INTEGER(first), *two = INTEGER(second);
    double *out = REAL(result);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) num_threads(threads) if(threads > 1)
#endif
    for (int k = 0; k < rank; ++k) {
        std::vector<int> rows;
        for (int s = 0; s < samples; ++s) {
            const Sample &x = input[s];
            const size_t pos = k + static_cast<size_t>(rank) * s;
            const double *h = x.profile + static_cast<size_t>(k) * x.scans;
            rows.clear();
            for (int i = 0; i < x.scans; ++i) if (h[i] > 0) rows.push_back(i);
            if (rows.empty()) { out[pos] = NA_REAL; continue; }
            const int j = one[pos] - 1, l = two[pos] - 1;
            double remove = 0, single = 0, blend = 0;
            for (int m = 0; m < ions; ++m) {
                const double *column = p + static_cast<size_t>(m) * rank;
                const double a = -column[k];
                const double b = j >= 0 ? a + column[j] : 0;
                const double c = l >= 0 && j >= 0 ?
                    a + (1 - weights[pos]) * column[j] + weights[pos] * column[l] : 0;
                if (a == 0 && b == 0 && c == 0) continue;
                const double *y = x.observed + static_cast<size_t>(m) * x.scans;
                const double *mu = x.mean + static_cast<size_t>(m) * x.scans;
                for (int i : rows) {
                    if (a != 0) remove += loss_change(y[i], mu[i], a * h[i]);
                    if (b != 0) single += loss_change(y[i], mu[i], b * h[i]);
                    if (c != 0) blend += loss_change(y[i], mu[i], c * h[i]);
                }
            }
            double best = j >= 0 ? std::min(remove, single) : remove;
            if (l >= 0 && j >= 0) best = std::min(best, blend);
            out[pos] = 2 * best / x.dispersion;
        }
    }
    return result;
}
