#include <Rcpp.h>
#include "count_profile.h"
#include "spectral_coordinate.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using logcon_counts::Vec;
using logcon_counts::Window;
namespace {
struct View {
    double *ptr;
    int rows, cols;
    explicit View(NumericMatrix x) : ptr(REAL(x)), rows(x.nrow()), cols(x.ncol()) {}
    double &operator()(int i, int j) const { return ptr[i + (size_t)rows * j]; }
};
Window window(View y, View o) {
    if (y.rows != o.rows || y.cols != o.cols)
        stop("Incompatible original-ion window.");
    Window w(y.rows, y.cols);
    for (int i = 0; i < y.rows; ++i)
        for (int j = 0; j < y.cols; ++j) {
            size_t k = (size_t)i * y.cols + j;
            w.y[k] = y(i, j);
            w.o[k] = o(i, j);
        }
    return w;
}
Vec values(NumericVector v) {
    return Vec(v.begin(), v.end());
}
struct Sample {
    View y, bg, mean, h, z, alpha;
    Vec t, warp;
    double phi, width;
    explicit Sample(Environment s)
        : y(as<NumericMatrix>(s["y"])), bg(as<NumericMatrix>(s["bg"])),
          mean(as<NumericMatrix>(s["mu"])), h(as<NumericMatrix>(s["h"])),
          z(as<NumericMatrix>(s["z"])), alpha(as<NumericMatrix>(s["alpha"])), t(values(s["t"])),
          warp(values(s["warp"])), phi(as<double>(s["phi"])), width(as<double>(s["width"])) {}
};
std::vector<Sample> samples(List a) {
    std::vector<Sample> out;
    out.reserve(a.size());
    for (SEXP s : a)
        out.emplace_back(Environment(s));
    return out;
}
std::vector<int> indexes(IntegerVector x) {
    std::vector<int> v(x.begin(), x.end());
    for (int &k : v)
        --k;
    return v;
}
} // namespace
// [[Rcpp::export]]
NumericVector cc_log_shape(NumericVector time, NumericVector h) {
    if (time.size() != h.size())
        stop("Unequal time and profile lengths.");
    return wrap(logcon_counts::shape(values(time), values(h)));
}
// [[Rcpp::export]]
IntegerVector cc_faces(NumericVector time, NumericMatrix h, bool validate = false) {
    const Vec t = values(time);
    IntegerVector result(h.ncol());
    if ((int)t.size() != h.nrow())
        stop("Unequal time and profile lengths.");
    for (int j = 0; j < h.ncol(); ++j) {
        Vec x(REAL(h) + (size_t)j * h.nrow(), REAL(h) + (size_t)(j + 1) * h.nrow());
        double err;
        result[j] = logcon_counts::face(t, x, &err);
        if (validate && err > 1e-6)
            stop("Nonconcave fitted profile on native scan coordinates.");
    }
    return result;
}
// [[Rcpp::export]]
NumericVector cc_local_profile(NumericVector t, NumericMatrix y, NumericMatrix o, NumericVector p,
                               NumericVector h, double phi, double price, bool scale_only = false) {
    Window w = window(View(y), View(o));
    if (t.size() != w.n || h.size() != w.n || p.size() != w.m)
        stop("Incompatible profile window.");
    double delta = 0;
    return wrap(scale_only ? w.rescale(values(p), values(h))
                           : w.update(values(t), values(p), values(h), phi, price, delta));
}
// [[Rcpp::export]]
NumericVector cc_local_spectrum(NumericMatrix y, NumericMatrix o, NumericVector h,
                                NumericVector old, double price, int threads = 1) {
    const int n = y.nrow(), m = y.ncol();
    if (y.nrow() != o.nrow() || y.ncol() != o.ncol() || h.size() != n || old.size() != m ||
        threads < 1)
        stop("Incompatible spectrum window.");
    const double *yp = REAL(y), *op = REAL(o), *hp = REAL(h), *pr = REAL(old);
    NumericVector result(m);
    double *out = REAL(result), total = std::accumulate(h.begin(), h.end(), 0.);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)                                                          \
    num_threads(threads) if (threads > 1 && (double)n * m > 50000)
#endif
    for (int j = 0; j < m; ++j) {
        size_t i = (size_t)n * j;
        out[j] = logcon_counts::solve_ion(yp + i, op + i, hp, n, total, pr[j], price, false).value;
    }
    return result;
}
// [[Rcpp::export]]
double cc_profile_sweep(List data, NumericMatrix spectra, IntegerVector order, NumericMatrix bands,
                        double price, bool birth, int threads) {
    auto s = samples(data);
    View p(spectra), b(bands);
    auto ord = indexes(order);
    Vec change(s.size());
    const int ns = s.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threads) if (threads > 1 && ns > 1)
#endif
    for (int si = 0; si < ns; ++si) {
        auto &x = s[si];
        const int n = x.h.rows;
        for (int k : ord) {
            std::vector<int> ion;
            Vec pc;
            for (int j = 0; j < p.cols; ++j)
                if (p(k, j) > 0) {
                    ion.push_back(j);
                    pc.push_back(p(k, j));
                }
            if (ion.empty())
                continue;
            int lo = n, hi = 0;
            for (int i = 0; i < n; ++i)
                if (x.h(i, k) > 0) {
                    lo = std::min(lo, i);
                    hi = i + 1;
                }
            if (hi == 0 && !birth)
                continue;
            lo = std::min(lo,
                          (int)(std::lower_bound(x.warp.begin(), x.warp.end(), b(k, 0) - x.width) -
                                x.warp.begin()));
            hi = std::max(hi,
                          (int)(std::lower_bound(x.warp.begin(), x.warp.end(), b(k, 1) + x.width) -
                                x.warp.begin()));
            if (lo >= hi)
                continue;
            for (;;) {
                Window w(hi - lo, ion.size());
                Vec prev(w.n), t = logcon_counts::slice(x.t, lo, hi);
                for (int i = 0; i < w.n; ++i) {
                    prev[i] = x.h(lo + i, k);
                    for (int j = 0; j < w.m; ++j) {
                        const size_t ij = (size_t)i * w.m + j;
                        w.y[ij] = x.y(lo + i, ion[j]);
                        w.o[ij] = std::max(1e-8, x.mean(lo + i, ion[j]) - prev[i] * pc[j]);
                    }
                }
                double gain;
                Vec v = w.update(t, pc, prev, x.phi, price, gain);
                double cutoff = .01 * logcon_counts::maximum(v);
                int extend = std::max(3, (hi - lo) / 2);
                int l = v.front() > cutoff ? std::max(0, lo - extend) : lo,
                    r = v.back() > cutoff ? std::min(n, hi + extend) : hi;
                if (l != lo || r != hi) {
                    lo = l;
                    hi = r;
                    continue;
                }
                change[si] +=
                    gain + price * (logcon_counts::face(t, v) - logcon_counts::face(t, prev));
                for (int i = 0; i < w.n; ++i) {
                    x.h(lo + i, k) = v[i];
                    for (int j = 0; j < w.m; ++j)
                        x.mean(lo + i, ion[j]) = w.o[(size_t)i * w.m + j] + v[i] * pc[j];
                }
                break;
            }
        }
    }
    return logcon_counts::sum(change);
}
// A component keeps its Gauss-Seidel order within each ion. Different ion columns are independent.
// [[Rcpp::export]]
double cc_spectrum_sweep(List data, NumericMatrix spectra, IntegerVector order, double price,
                         bool track, int threads) {
    auto s = samples(data);
    View p(spectra);
    auto ord = indexes(order);
    struct Entry {
        int si, i;
        double raw, scaled;
    };
    std::vector<std::vector<Entry>> positions(p.rows);
    for (int k : ord)
        for (int si = 0; si < (int)s.size(); ++si)
            for (int i = 0; i < s[si].h.rows; ++i) {
                double h = s[si].h(i, k);
                if (h > 0)
                    positions[k].push_back({si, i, h, h / s[si].phi});
            }
    Vec gains(p.cols);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads) if (threads > 1)
#endif
    for (int j = 0; j < p.cols; ++j) {
        Vec y, o, h, raw;
        for (int k : ord) {
            auto &at = positions[k];
            int n = at.size();
            if (!n)
                continue;
            y.resize(n);
            o.resize(n);
            h.resize(n);
            raw.resize(n);
            double old = p(k, j);
            for (int i = 0; i < n; ++i) {
                auto &e = at[i];
                auto &x = s[e.si];
                h[i] = e.scaled;
                y[i] = x.y(e.i, j) / x.phi;
                raw[i] = std::max(1e-8, x.mean(e.i, j) - e.raw * old);
                o[i] = raw[i] / x.phi;
            }
            auto fit = logcon_counts::solve_ion(y.data(), o.data(), h.data(), n,
                                                logcon_counts::sum(h), old, price, track);
            p(k, j) = fit.value;
            gains[j] += fit.gain;
            for (int i = 0; i < n; ++i) {
                auto &e = at[i];
                s[e.si].mean(e.i, j) = raw[i] + e.raw * fit.value;
            }
        }
        // Background coefficients are unpenalised, sample-specific nonnegative fits.
        for (auto &x : s)
            for (int k = 0; k < x.z.cols; ++k) {
                std::vector<int> at;
                y.clear();
                o.clear();
                h.clear();
                double old = x.alpha(k, j);
                for (int i = 0; i < x.z.rows; ++i)
                    if (x.z(i, k) > 0) {
                        at.push_back(i);
                        y.push_back(x.y(i, j) / x.phi);
                        h.push_back(x.z(i, k) / x.phi);
                        o.push_back(std::max(1e-8, x.mean(i, j) - x.z(i, k) * old) / x.phi);
                    }
                auto f = logcon_counts::solve_ion(y.data(), o.data(), h.data(), at.size(),
                                                  logcon_counts::sum(h), old, 0., track);
                for (int i : at) {
                    double dx = x.z(i, k) * (f.value - old);
                    x.mean(i, j) += dx;
                    x.bg(i, j) += dx;
                }
                x.alpha(k, j) = f.value;
                gains[j] += f.gain;
            }
        for (auto &x : s)
            for (int i = 0; i < x.y.rows; ++i) {
                x.bg(i, j) = std::max(x.bg(i, j), 1e-8);
                x.mean(i, j) = std::max(x.mean(i, j), 1e-8);
            }
    }
    for (int k : ord)
        if (!positions[k].empty()) {
            double total = 0;
            for (int j = 0; j < p.cols; ++j)
                total += p(k, j);
            for (auto &e : positions[k])
                s[e.si].h(e.i, k) *= total;
            if (total > 0)
                for (int j = 0; j < p.cols; ++j)
                    p(k, j) /= total;
        }
    return 2 * logcon_counts::sum(gains);
}
// [[Rcpp::export]]
double cc_deviance(NumericMatrix y, NumericMatrix mu, double phi) {
    if (y.nrow() != mu.nrow() || y.ncol() != mu.ncol() || !(phi > 0))
        stop("Invalid count deviance dimensions or dispersion.");
    double total = 0;
    for (R_xlen_t i = 0; i < y.size(); ++i) {
        double v = y[i], m = mu[i];
        total += (v > 0 ? v * std::log(v / m) : 0) - v + m;
    }
    return 2 * total / phi;
}
namespace {
inline double count_loss_change(double y, double mu, double old) {
    const double change = mu - old, relative = change / old;
    const double log_ratio =
        std::abs(relative) < 0.5 ? std::log1p(relative) : std::log(mu) - std::log(old);
    return change - y * log_ratio;
}
} // namespace
// Change in original-channel Poisson or quasi-Poisson deviance on a local window.
// Compute the log ratio stably without materializing a ratio or logarithm matrix.
// [[Rcpp::export]]
double cc_diff(NumericMatrix y, NumericMatrix mu, NumericMatrix old, double phi) {
    if (y.nrow() != mu.nrow() || y.ncol() != mu.ncol() || old.nrow() != mu.nrow() ||
        old.ncol() != mu.ncol() || !(phi > 0))
        stop("Incompatible local count-deviance matrices or dispersion.");
    long double total = 0;
    for (R_xlen_t i = 0; i < y.size(); ++i)
        total += count_loss_change(y[i], mu[i], old[i]);
    return (2 / phi) * static_cast<double>(total);
}
// Extract the enclosing native scan interval for one local factor update.
// Only its observed counts, focal profiles and offsets are materialized.
// [[Rcpp::export]]
SEXP cc_local_window(Environment sample, NumericMatrix spectra, IntegerVector components,
                     IntegerVector ions, int si, bool floor_background) {
    const View h(as<NumericMatrix>(sample["h"])), p(spectra);
    int left = h.rows, right = -1;
    for (int column : components) {
        if (column < 1 || column > h.cols || column > p.rows)
            stop("Component index outside the local fit.");
        const double *v = h.ptr + static_cast<size_t>(column - 1) * h.rows;
        for (int i = 0; i < left; ++i)
            if (v[i] > 0) {
                left = i;
                break;
            }
        for (int i = h.rows - 1; i > right; --i)
            if (v[i] > 0) {
                right = i;
                break;
            }
    }
    if (right < left)
        return R_NilValue;
    const int n = right - left + 1, m = ions.size(), k = components.size();
    const View y(as<NumericMatrix>(sample["y"])), mean(as<NumericMatrix>(sample["mu"])),
        background(as<NumericMatrix>(sample["bg"]));
    const NumericVector time = sample["t"];
    if (time.size() != h.rows || y.rows != h.rows || mean.rows != y.rows ||
        background.rows != y.rows || mean.cols != y.cols || background.cols != y.cols ||
        p.cols != y.cols)
        stop("Incompatible sample matrices for a local fit.");
    NumericMatrix observed(n, m), old(n, m), off(n, m), profiles(n, k);
    for (int j = 0; j < k; ++j)
        std::copy_n(h.ptr + static_cast<size_t>(components[j] - 1) * h.rows + left, n,
                    profiles.begin() + static_cast<size_t>(j) * n);
    for (int j = 0; j < m; ++j) {
        const int ion = ions[j] - 1;
        if (ion < 0 || ion >= y.cols)
            stop("Ion index outside the local fit.");
        const size_t offset = static_cast<size_t>(ion) * y.rows + left,
                     dest = static_cast<size_t>(j) * n;
        std::copy_n(y.ptr + offset, n, observed.begin() + dest);
        std::copy_n(mean.ptr + offset, n, old.begin() + dest);
        for (int i = 0; i < n; ++i) {
            double focal = 0;
            for (int g = 0; g < k; ++g)
                focal += profiles(i, g) * p(components[g] - 1, ion);
            off(i, j) =
                std::max(old(i, j) - focal, floor_background ? background(left + i, ion) : 1e-8);
        }
    }
    return List::create(_["si"] = si, _["r"] = seq(left + 1, right + 1),
                        _["t"] = NumericVector(time.begin() + left, time.begin() + right + 1),
                        _["y"] = observed, _["off"] = off, _["old"] = old, _["h"] = profiles,
                        _["phi"] = sample["phi"]);
}
// Remove sample-level profiles or an entire shared spectrum when that lowers the
// same original-ion criterion. The caller orders components by their current area.
// [[Rcpp::export]]
double cc_support_sweep(List data, NumericMatrix spectra, IntegerVector order, double price) {
    auto s = samples(data);
    View p(spectra);
    double gain = 0;
    for (int k : indexes(order)) {
        std::vector<int> ions;
        if (k < 0 || k >= p.rows)
            stop("Component index outside the support update.");
        for (int j = 0; j < p.cols; ++j)
            if (p(k, j) > 0)
                ions.push_back(j);
        if (ions.empty())
            continue;
        std::vector<std::vector<int>> active(s.size());
        Vec dev(s.size()), merit(s.size());
        long double benefit = 0;
        for (size_t i = 0; i < s.size(); ++i) {
            auto &x = s[i];
            Vec t, h;
            for (int row = 0; row < x.h.rows; ++row)
                if (x.h(row, k) > 0) {
                    active[i].push_back(row);
                    t.push_back(x.t[row]);
                    h.push_back(x.h(row, k));
                }
            if (h.empty())
                continue;
            long double delta = 0;
            for (int j : ions)
                for (int row : active[i]) {
                    double m = x.mean(row, j);
                    double o = std::max(1e-8, m - x.h(row, k) * p(k, j));
                    delta += count_loss_change(x.y(row, j), o, m);
                }
            dev[i] = (2 / x.phi) * static_cast<double>(delta);
            merit[i] = dev[i] - price * logcon_counts::face(t, h);
            benefit += std::max(merit[i], 0.);
        }
        bool retire = static_cast<double>(benefit) <= price * (ions.size() - 1);
        for (size_t i = 0; i < s.size(); ++i)
            if (retire || merit[i] < 0) {
                auto &x = s[i];
                for (int row : active[i]) {
                    double h = x.h(row, k);
                    for (int j : ions)
                        x.mean(row, j) = std::max(1e-8, x.mean(row, j) - h * p(k, j));
                    x.h(row, k) = 0;
                }
                gain -= dev[i];
            }
        if (retire)
            for (int j : ions)
                p(k, j) = 0;
    }
    return gain;
}
// Mutate only the fitter's owned matrices; local trials keep copies until accepted.
// [[Rcpp::export]]
void cc_write_block(NumericMatrix dest, IntegerVector rows, IntegerVector columns,
                    NumericMatrix block) {
    if (rows.size() != block.nrow() || columns.size() != block.ncol())
        stop("Invalid update dimensions.");
    for (int i : rows)
        if (i < 1 || i > dest.nrow())
            stop("Scan index outside the local update.");
    for (int j : columns)
        if (j < 1 || j > dest.ncol())
            stop("Channel or component index outside the local update.");
    for (int j = 0; j < block.ncol(); ++j)
        for (int i = 0; i < block.nrow(); ++i)
            dest(rows[i] - 1, columns[j] - 1) = block(i, j);
}
