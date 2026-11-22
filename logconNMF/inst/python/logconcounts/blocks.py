"""Coupled count updates on disjoint samples or ion columns."""

import numpy as np
from numba import njit
from scipy import sparse

from .kernels import face, solve_ion, update_h
from .workers import capacity, ordered_map


@njit(cache=True, nogil=True)
def profile_sample(t, warp, y, phi, h, mu, P, order, bands, width, price, birth):
    gains = np.zeros(len(P))
    for k in order:
        ion = np.flatnonzero(P[k] > 0)
        if not len(ion):
            continue
        old = h[:, k]
        present = np.flatnonzero(old > 0)
        if not len(present) and not birth:
            continue
        a = np.searchsorted(warp, bands[k, 0] - width)
        b = np.searchsorted(warp, bands[k, 1] + width)
        if len(present):
            a = min(a, present[0])
            b = max(b, present[-1] + 1)
        if b <= a:
            continue
        pc = P[k, ion]
        while True:
            prev = old[a:b].copy()
            offsets = np.empty((b - a, len(ion)))
            counts = np.empty_like(offsets)
            for row in range(b - a):
                for j, m in enumerate(ion):
                    counts[row, j] = y[a + row, m]
                    offsets[row, j] = max(1e-8, mu[a + row, m] - prev[row] * pc[j])
            v, gain = update_h(t[a:b], counts, offsets, pc, prev, phi, price)
            peak = v.max()
            extension = max(3, (b - a) // 2)
            lower = max(0, a - extension) if v[0] > 0.01 * peak else a
            upper = min(len(t), b + extension) if v[-1] > 0.01 * peak else b
            if lower == a and upper == b:
                break
            a, b = lower, upper
        gains[k] = gain + price * (face(t[a:b], v) - face(t[a:b], prev))
        old[a:b] = v
        for row in range(b - a):
            for j, m in enumerate(ion):
                mu[a + row, m] = offsets[row, j] + v[row] * pc[j]
    return gains


@njit(cache=True, nogil=True)
def ion_sweep(
    lo, hi, order, ptr, sample, row, raw_h, h, phi, y, mu, P, q, gain, price, track
):
    longest = max(1, np.max(np.diff(ptr)))
    counts = np.empty(longest)
    offsets = np.empty(longest)
    unscaled = np.empty(longest)
    for k in order:
        a, b = ptr[k], ptr[k + 1]
        n = b - a
        if n == 0:
            continue
        hv = h[a:b]
        total = hv.sum()
        for m in range(lo, hi):
            old = P[k, m]
            for j in range(n):
                i = a + j
                s, r = sample[i], row[i]
                counts[j] = y[s][r, m] / phi[s]
                unscaled[j] = max(1e-8, mu[s][r, m] - raw_h[i] * old)
                offsets[j] = unscaled[j] / phi[s]
            coef, delta = solve_ion(
                counts[:n], offsets[:n], hv, total, old, price, track
            )
            q[k, m] = coef
            gain[k, m] = delta
            for j in range(n):
                i = a + j
                mu[sample[i]][row[i], m] = unscaled[j] + raw_h[i] * coef


def update_spectra(data, H, P, mu, price, track):
    # Within an ion, preserve the original Gauss-Seidel order. Its other ions
    # enter only in the unit-sum rescaling, which can wait until the entire sweep.
    mass = np.sum([h.sum(0) for h in H], axis=0)
    if not np.any(mass > 0):
        return 0.0
    order = np.argsort(-mass)
    design = sparse.vstack([sparse.csc_matrix(h) for h in H], format="csc")
    ends = np.cumsum([len(h) for h in H])
    sample = np.searchsorted(ends, design.indices, side="right")
    row = design.indices - np.r_[0, ends[:-1]][sample]
    phi = np.array([s.phi for s in data])
    hs = design.data / phi[sample]
    # Numba restricts homogeneous tuple arguments to at most 1,000 entries.
    if len(data) > 1000:
        from numba.typed import List

        y = List([s.y for s in data])
        mus = List(mu)
    else:
        y = tuple(s.y for s in data)
        mus = tuple(mu)
    q = P.copy()
    gains = np.zeros_like(q)
    width = max(1, (P.shape[1] + capacity() - 1) // capacity())

    def sweep(j):
        ion_sweep(
            j,
            min(j + width, P.shape[1]),
            order,
            design.indptr,
            sample,
            row,
            design.data,
            hs,
            phi,
            y,
            mus,
            P,
            q,
            gains,
            price,
            track,
        )

    ordered_map(sweep, range(0, P.shape[1], width))
    total = q.sum(axis=1)
    live = mass > 0
    for h in H:
        h *= np.where(live, total, 1.0)
    P[live] = q[live] / np.where(total[live] > 0, total[live], 1.0)[:, None]
    return 2 * gains.sum()
