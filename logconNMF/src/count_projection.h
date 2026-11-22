#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
namespace logcon_geometry {
using Vec = std::vector<double>;
inline double largest(const Vec &x) {
    return *std::max_element(x.begin(), x.end());
}

inline Vec fit_knots(const Vec &time, const Vec &y, const Vec &w, const std::vector<int> &knots) {
    const int n = time.size(), p = knots.size();
    Vec diagonal(p), upper(p), rhs(p), z(n);
    for (int i = 0, k = 0; i < n; ++i) {
        if (k + 2 < p && i >= knots[k + 1])
            ++k;
        double b = (time[i] - time[knots[k]]) / (time[knots[k + 1]] - time[knots[k]]), a = 1 - b;
        diagonal[k] += w[i] * a * a;
        diagonal[k + 1] += w[i] * b * b;
        upper[k] += w[i] * a * b;
        rhs[k] += w[i] * a * y[i];
        rhs[k + 1] += w[i] * b * y[i];
    }
    for (int k = 1; k < p; ++k) {
        const double r = upper[k - 1] / diagonal[k - 1];
        diagonal[k] -= r * upper[k - 1];
        rhs[k] -= r * rhs[k - 1];
    }
    rhs.back() /= diagonal.back();
    for (int k = p - 2; k >= 0; --k)
        rhs[k] = (rhs[k] - upper[k] * rhs[k + 1]) / diagonal[k];
    for (int i = 0, k = 0; i < n; ++i) {
        if (k + 2 < p && i >= knots[k + 1])
            ++k;
        const double b = (time[i] - time[knots[k]]) / (time[knots[k + 1]] - time[knots[k]]);
        z[i] = (1 - b) * rhs[k] + b * rhs[k + 1];
    }
    return z;
}

inline Vec project(const Vec &t, const Vec &y, const Vec &w) {
    const int n = t.size();
    if (n < 3)
        return y;
    std::vector<int> knots{0, n - 1};
    Vec z = fit_knots(t, y, w, knots);
    for (int it = 0; it < 10 * n + 10; ++it) {
        Vec candidate = fit_knots(t, y, w, knots);
        int drop = -1;
        double fraction = 1;
        for (int j = 1; j + 1 < (int)knots.size(); ++j) {
            const int l = knots[j - 1], m = knots[j], r = knots[j + 1];
            const double b = (t[m] - t[l]) / (t[r] - t[l]);
            double slack = candidate[m] - (1 - b) * candidate[l] - b * candidate[r];
            if (slack < -1e-10) {
                double old = std::max(0., z[m] - (1 - b) * z[l] - b * z[r]);
                double step = old / (old - slack);
                if (step < fraction) {
                    fraction = step;
                    drop = j;
                }
            }
        }
        if (drop >= 0) {
            for (int i = 0; i < n; ++i)
                z[i] += fraction * (candidate[i] - z[i]);
            knots.erase(knots.begin() + drop);
            continue;
        }
        z.swap(candidate);
        double total = 0, moment = 0, norm = 0;
        for (int i = 0; i < n; ++i)
            norm += std::abs(w[i] * (z[i] - y[i]));
        double best = 1e-12 * std::max(1., norm);
        int add = -1;
        // Directional derivative for inserting each inactive hinge.
        for (int i = n - 1; i > 0; --i) {
            const double g = w[i] * (z[i] - y[i]);
            total += g;
            moment += t[i] * g;
            double score = moment - t[i - 1] * total;
            if (score > best && !std::binary_search(knots.begin(), knots.end(), i - 1)) {
                best = score;
                add = i - 1;
            }
        }
        if (add < 0)
            return z;
        knots.insert(std::lower_bound(knots.begin(), knots.end(), add), add);
    }
    // Last feasible iterate is a safe proposal; the likelihood decides acceptance.
    return z;
}

} // namespace logcon_geometry
