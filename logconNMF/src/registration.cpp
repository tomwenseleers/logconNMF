#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;

// Find the 1%-height support and its interpolated onset directly on the native grid.
// The dictionary window includes exact threshold ties; the alignment width uses the
// open level set, as in the shape definition. O(nnz) for a compressed sparse matrix.
// [[Rcpp::export]]
List slcn_profile_landmarks(SEXP H, NumericVector time) {
    const bool sparse = Rf_isS4(H);
    IntegerVector ptr, rid;
    NumericVector values;
    int rows, k;
    if (sparse) {
        S4 h(H);
        IntegerVector dim = h.slot("Dim");
        rows = dim[0];
        k = dim[1];
        ptr = h.slot("p");
        rid = h.slot("i");
        values = h.slot("x");
    } else {
        NumericMatrix h = as<NumericMatrix>(H);
        rows = h.nrow();
        k = h.ncol();
        values = NumericVector(static_cast<SEXP>(h));
    }
    if (time.size() != rows)
        stop("Profile row count differs from the time grid.");
    IntegerMatrix windows(k, 3);
    NumericVector start(k), apex(k), width(k);
    for (int j = 0; j < k; ++j) {
        const R_xlen_t begin = sparse ? ptr[j] : static_cast<R_xlen_t>(j) * rows;
        const R_xlen_t end = sparse ? ptr[j + 1] : begin + rows;
        if (begin == end)
            stop("A profile has empty support.");
        const double *v = REAL(values) + begin;
        const int len = static_cast<int>(end - begin);
        const int peak = std::max_element(v, v + len) - v;
        const double top = v[peak], threshold = 0.01 * top;
        if (!(top > 0))
            stop("A profile has empty support.");
        auto row = [&](int i) { return sparse ? rid[begin + i] : i; };
        int lo = 0, hi = len - 1;
        while (v[lo] < threshold)
            ++lo;
        while (v[hi] < threshold)
            --hi;
        windows(j, 0) = row(lo) + 1;
        windows(j, 1) = row(hi) + 1;
        windows(j, 2) = row(peak) + 1;
        while (v[lo] <= threshold)
            ++lo;
        while (v[hi] <= threshold)
            --hi;
        const int first = row(lo);
        double on = time[first];
        if (first > 0) {
            const double prev = lo > 0 && row(lo - 1) == first - 1 ? v[lo - 1] : 0.;
            on = time[first - 1] +
                 (time[first] - time[first - 1]) * (threshold - prev) / (v[lo] - prev);
        }
        start[j] = on;
        apex[j] = time[row(peak)];
        width[j] = time[row(hi)] - time[first];
    }
    return List::create(_["windows"] = windows, _["start"] = start, _["apex"] = apex,
                        _["width"] = width);
}

namespace {
struct Axis {
    std::vector<double> stored, norms;
    int rows, cols;
    Axis(NumericMatrix s, bool hellinger)
        : stored(s.size()), norms(s.nrow()), rows(s.nrow()), cols(s.ncol()) {
        // R's rowSums accumulates in extended precision where the platform supports it.
        for (int k = 0; k < rows; ++k) {
            long double sum = 0.;
            for (int j = 0; j < cols; ++j) {
                double v = s(k, j);
                sum += hellinger ? v : v * v;
            }
            norms[k] = static_cast<double>(sum);
            for (int j = 0; j < cols; ++j)
                stored[static_cast<size_t>(k) * cols + j] =
                    hellinger ? std::sqrt(s(k, j) / norms[k]) : s(k, j);
        }
    }
    const double *row(int i) const { return stored.data() + static_cast<size_t>(i) * cols; }
};
} // namespace

// Nearest spectral neighbour in both directions, restricted to eligible retention-time pairs.
// Return reciprocal pairs in the original sample order. Ties follow the original row order.
// [[Rcpp::export]]
List slcn_anchor_matches(NumericMatrix s, NumericMatrix ref, NumericVector t, NumericVector rt,
                         double window, bool hellinger) {
    if (s.ncol() != ref.ncol() || s.nrow() != t.size() || ref.nrow() != rt.size() || window < 0)
        stop("Incompatible spectral and retention-time axes.");
    Axis a(s, hellinger), b(ref, hellinger);
    std::vector<int> sorted;
    for (int i = 0; i < b.rows; ++i)
        if (R_FINITE(rt[i]))
            sorted.push_back(i);
    std::stable_sort(sorted.begin(), sorted.end(), [&](int i, int j) { return rt[i] < rt[j]; });
    std::vector<int> forward(a.rows, -1), reverse(b.rows, -1);
    std::vector<double> fd(a.rows, R_PosInf), rd(b.rows, R_PosInf);
    for (int i = 0; i < a.rows; ++i) {
        if (!R_FINITE(t[i]))
            continue;
        double lower = std::nextafter(t[i] - window, -std::numeric_limits<double>::infinity());
        double upper = std::nextafter(t[i] + window, std::numeric_limits<double>::infinity());
        auto first = std::lower_bound(sorted.begin(), sorted.end(), lower,
                                      [&](int j, double x) { return rt[j] < x; });
        const double *x = a.row(i);
        for (auto at = first; at != sorted.end() && rt[*at] <= upper; ++at) {
            int j = *at;
            if (std::abs(t[i] - rt[j]) > window)
                continue;
            double prod = std::inner_product(x, x + a.cols, b.row(j), 0.);
            if (!R_FINITE(prod))
                continue;
            double d = hellinger ? std::sqrt(std::max(0., 1. - prod))
                                 : 1. - prod / std::sqrt(a.norms[i] * b.norms[j]);
            if (d < fd[i] || (d == fd[i] && j < forward[i])) {
                fd[i] = d;
                forward[i] = j;
            }
            if (d < rd[j]) {
                rd[j] = d;
                reverse[j] = i;
            }
        }
    }
    std::vector<int> from, to;
    std::vector<double> distances;
    for (int i = 0; i < a.rows; ++i)
        if (forward[i] >= 0 && reverse[forward[i]] == i) {
            from.push_back(i + 1);
            to.push_back(forward[i] + 1);
            distances.push_back(fd[i]);
        }
    return List::create(_["sample"] = from, _["reference"] = to, _["distance"] = distances);
}
