#pragma once
#include "count_projection.h"
#include <limits>

namespace logcon_counts {
using logcon_geometry::Vec;
inline Vec slice(const Vec &v, int first, int end) {
    return Vec(v.begin() + first, v.begin() + end);
}
inline double sum(const Vec &v) {
    return std::accumulate(v.begin(), v.end(), 0.);
}
inline double maximum(const Vec &v) {
    return v.empty() ? 0. : *std::max_element(v.begin(), v.end());
}

// All shape constraints use physical scan coordinates, including on irregular grids.
inline Vec project_log(Vec time, Vec y, Vec w) {
    if (time.size() <= 2)
        return y;
    const double start = time.front(), width = time.back() - start, scale = maximum(w);
    for (int i = 0; i < (int)y.size(); ++i) {
        time[i] = (time[i] - start) / width;
        w[i] /= scale;
    }
    return logcon_geometry::project(time, y, w);
}
inline Vec shape(const Vec &t, const Vec &h) {
    const int n = h.size();
    int lo = 0, hi = n;
    while (lo < hi && h[lo] <= 0)
        ++lo;
    while (lo < hi && h[hi - 1] <= 0)
        --hi;
    Vec result(n);
    if (lo == hi)
        return result;
    const double peak = maximum(h);
    Vec y(hi - lo), w(hi - lo);
    for (int i = lo; i < hi; ++i) {
        y[i - lo] = std::log(std::max(h[i], peak * 1e-9));
        w[i - lo] = std::max(h[i] / peak, 1e-6);
    }
    Vec z = project_log(slice(t, lo, hi), y, w);
    for (int i = lo; i < hi; ++i)
        result[i] = std::exp(std::min(z[i - lo], std::log(peak) + 4.));
    return result;
}
inline int face(const Vec &t, const Vec &h, double *violation = nullptr) {
    int n = h.size(), count = 0, prev = -1, last = -1, df = 0;
    double maxlog = 0, error = 0;
    for (double v : h)
        if (v > 0)
            maxlog = std::max(maxlog, std::abs(std::log(v)));
    const double tol = std::sqrt(std::numeric_limits<double>::epsilon()) * (1 + maxlog);
    for (int i = 0; i < n; ++i)
        if (h[i] > 0) {
            if (last >= 0 && i != last + 1)
                error = std::numeric_limits<double>::infinity();
            if (prev >= 0) {
                double left = std::log(h[prev]);
                double gap = std::log(h[last]) - left -
                             (std::log(h[i]) - left) * (t[last] - t[prev]) / (t[i] - t[prev]);
                if (gap > tol)
                    ++df;
                error = std::max(error, -gap);
            }
            prev = last;
            last = i;
            ++count;
        }
    if (violation)
        *violation = error;
    return std::min(count, 2) + df;
}

// Local matrices in scan order, contiguous over ions.
struct Window {
    int n, m;
    Vec y, o;
    Window(int rows, int ions) : n(rows), m(ions), y((size_t)rows * ions), o(y.size()) {}
    double gain(const Vec &h, const Vec &p, double phi) const {
        double g = 0;
        for (int i = 0; i < n; ++i)
            if (h[i] > 0)
                for (int j = 0; j < m; ++j)
                    if (p[j] > 0) {
                        size_t k = (size_t)i * m + j;
                        double v = h[i] * p[j];
                        g += y[k] * std::log1p(v / o[k]) - v;
                    }
        return 2 * g / phi;
    }
    Vec rescale(const Vec &p, const Vec &h) const {
        const double hp = maximum(h);
        if (hp <= 0)
            return h;
        Vec v = h;
        for (double &a : v)
            a /= hp;
        const double norm = sum(p) * sum(v);
        if (norm <= 0)
            return Vec(n);
        double lo = 0, hi = sum(y) / norm, a = std::min(hi, hp);
        for (int it = 0; it < 18; ++it) {
            double g = -norm, c = 0;
            for (int i = 0; i < n; ++i)
                if (v[i] > 0)
                    for (int j = 0; j < m; ++j)
                        if (p[j] > 0) {
                            size_t k = (size_t)i * m + j;
                            double z = v[i] * p[j];
                            z /= o[k] + a * z;
                            g += y[k] * z;
                            c += y[k] * z * z;
                        }
            if (std::abs(g) <= 1e-8 * norm)
                break;
            if (g > 0)
                lo = a;
            else
                hi = a;
            if (hi <= hp * 1e-12)
                return Vec(n);
            double trial = a + g / std::max(c, 1e-60);
            a = (lo < trial && trial < hi) ? trial : 0.5 * (lo + hi);
        }
        for (double &x : v)
            x *= a;
        return v;
    }
    Vec update(const Vec &t, const Vec &p, const Vec &h, double phi, double price,
               double &gain_change) const {
        const double old = gain(h, p, phi) - price * face(t, h), sp = sum(p);
        double best_score = std::max(0., old);
        Vec best = old < 0 ? Vec(n) : h, v(n), info(n);
        int left = 0, right = 0, start = 0;
        double value = 0, cur = 0;
        for (int i = 0; i < n; ++i) {
            double g0 = -sp, total = 0, gscan = 0;
            const size_t off = (size_t)i * m;
            for (int j = 0; j < m; ++j) {
                g0 += p[j] * y[off + j] / o[off + j];
                total += y[off + j];
            }
            if (g0 > 0 && sp > 0) {
                double lo = 0, hi = total / sp, a = std::min(hi, std::max(h[i], .001 * hi));
                for (int it = 0; it < 19; ++it) {
                    double g = -sp, c = 0;
                    for (int j = 0; j < m; ++j) {
                        double z = p[j] / (o[off + j] + a * p[j]);
                        g += y[off + j] * z;
                        c += y[off + j] * z * z;
                    }
                    if (std::abs(g) < 1e-7 * sp)
                        break;
                    if (g > 0)
                        lo = a;
                    else
                        hi = a;
                    double trial = a + g / std::max(c, 1e-60);
                    a = (lo < trial && trial < hi) ? trial : .5 * (lo + hi);
                }
                v[i] = a;
                for (int j = 0; j < m; ++j) {
                    double term = a * p[j];
                    gscan += y[off + j] * std::log1p(term / o[off + j]) - term;
                    info[i] += p[j] * p[j] / (o[off + j] + term);
                }
            }
            cur += 2 * gscan / phi - price;
            if (cur > value) {
                value = cur;
                left = start;
                right = i + 1;
            }
            if (cur <= 0) {
                start = i + 1;
                cur = 0;
            }
        }
        auto consider = [&](const Vec &candidate) {
            Vec fit = rescale(p, candidate);
            double score = gain(fit, p, phi) - price * face(t, fit);
            if (score > best_score + 1e-7) {
                best_score = score;
                best.swap(fit);
            }
        };
        if (right > left) {
            left = std::max(0, left - 1);
            right = std::min(n, right + 1);
            Vec ylog = slice(v, left, right), w(right - left);
            double mx = maximum(ylog);
            if (mx > 0) {
                for (int j = 0; j < right - left; ++j) {
                    w[j] = std::max(info[left + j] * ylog[j] * ylog[j], 1e-15);
                    ylog[j] = std::log(std::max(ylog[j], mx * 1e-8));
                }
                double floor = maximum(w) * 1e-7;
                for (double &x : w)
                    x = std::max(x, floor);
                Vec z = project_log(slice(t, left, right), ylog, w), candidate(n);
                for (int i = left; i < right; ++i)
                    candidate[i] = std::exp(std::min(z[i - left], std::log(mx) + 3));
                consider(candidate);
            }
        }
        if (maximum(h) > 0)
            consider(h);
        gain_change = best_score - old;
        return best;
    }
};
} // namespace logcon_counts
