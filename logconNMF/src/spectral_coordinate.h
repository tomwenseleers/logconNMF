#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

// One ion's conditional count likelihood. Counts, offsets and profile have
// already been divided by the dispersion of their originating sample.
namespace logcon_counts {
struct Coordinate {
    double value;
    double gain;
};

inline Coordinate solve_ion(const double *counts, const double *offset, const double *h, int n,
                            double total, double previous, double price, bool track) {
    if (!(total > 0))
        return {0, 0};
    double at_zero = -total, count_sum = 0, previous_gain = 0;
    for (int i = 0; i < n; ++i) {
        count_sum += counts[i];
        at_zero += counts[i] * h[i] / offset[i];
        if (track && previous > 0)
            previous_gain += counts[i] * std::log1p(previous * h[i] / offset[i]);
    }
    previous_gain = track ? previous_gain - previous * total : 0;
    if (at_zero <= 0)
        return {0, -previous_gain};
    if (price > 0) {
        // For N counts and T=sum(h), Jensen bounds twice the maximum
        // log-likelihood gain by N*r^2/(1+r), r=score(0)/T. A coordinate
        // below its price can be set to zero without solving the score equation.
        const double r = at_zero / total;
        const double upper_r = r + 8. * (n + 1.) * std::numeric_limits<double>::epsilon() * (1. + r);
        if (count_sum * upper_r * (upper_r / (1. + upper_r)) < price)
            return {0, -previous_gain};
    }
    double lower = 0, upper = count_sum / total;
    double x = std::min(upper, std::max(previous, 1e-4 * upper));
    for (int step = 0; step < 25; ++step) {
        double score = -total, curvature = 0;
        for (int i = 0; i < n; ++i) {
            double r = h[i] / (offset[i] + x * h[i]);
            score += counts[i] * r;
            curvature += counts[i] * r * r;
        }
        if (std::abs(score) < 1e-8 * total)
            break;
        if (score > 0)
            lower = x;
        else
            upper = x;
        double next = x + score / std::max(curvature, 1e-90);
        x = (lower < next && next < upper) ? next : 0.5 * (lower + upper);
    }
    double gain = -x * total;
    for (int i = 0; i < n; ++i)
        gain += counts[i] * std::log1p(x * h[i] / offset[i]);
    if (2 * gain <= price)
        return {0, -previous_gain};
    return {x, track ? gain - previous_gain : 0};
}
} // namespace logcon_counts
