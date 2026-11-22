"""Local cone coding of PSF coefficients, followed by joint scan regression."""

import numpy as np
from numba import njit
from scipy.optimize import nnls
from .kernels import log_shape, update_p


@njit(cache=True)
def scan_unmix(Y, B, P, h, eligible, phi):
    for t in range(len(Y)):
        ix = np.flatnonzero(eligible[t])
        if not len(ix):
            continue
        b = h[t, ix].copy()
        p = P[ix]
        y = Y[t]
        for it in range(2):
            pred = np.maximum(b @ p + B[t], 0.0)
            w = 1 / np.maximum(y if it == 0 else pred, max(phi, 1.0))
            pw = p * w
            gram = pw @ p.T
            score = pw @ (y - B[t]) - gram @ b
            for sweep in range(32):
                change = 0.0
                for i in range(len(ix)):
                    a = max(0.0, b[i] + score[i] / max(gram[i, i], 1e-60)) - b[i]
                    b[i] += a
                    score -= gram[:, i] * a
                    change = max(change, abs(a))
                if change < 1e-6 * max(1.0, np.max(b)):
                    break
        h[t, ix] = b
    return h


def local_directions(data, A, price):
    N = len(A.p)
    raw = np.zeros_like(A.p)
    purity = np.zeros(N)
    for i in range(N):
        s = data[A.sample[i]]
        left, right, ap = s.grid[A.index[i]]
        core = max(1, round((right - left) / 6))
        sl = slice(max(0, ap - core), min(len(s.t), ap + core + 1))
        y = s.y[sl]
        o = s.bg[sl]
        r = np.maximum(y - o, 0.0)
        h = r.sum(1)
        z = update_p(y / s.phi, o / s.phi, h / s.phi, A.p[i], price)
        p = z / z.sum() if z.sum() > 0 else A.p[i]
        raw[i] = p
        # Rank-one residual energy is independent of units. The extra noise
        # term caps the priority of weak, baseline-dominated windows.
        d = r - h[:, None] * p
        err = (d * d).sum() + s.phi * y.sum()
        purity[i] = (r * r).sum() / max(err, 1.0)
    return raw, purity


def proposed_spectra(data, A, price, partial):
    if partial == "only":
        return [(i, A.p[i]) for i in np.argsort(-A.height)]
    raw, purity = local_directions(data, A, price)
    order = [(i, raw[i]) for i in np.argsort(-purity)]
    if partial:
        order.extend((i, A.p[i]) for i in np.argsort(-A.height))
    return order


def anchor_dictionary(data, A, price, partial=False):
    q = A.p.shape[1]
    P = []
    order = proposed_spectra(data, A, price, partial)
    bands = np.empty((len(order), 2))
    directions = []
    for i, p in order:
        u = p / max(np.linalg.norm(p), 1e-60)
        if len(P):
            local = bands[:len(P)]
            c = np.flatnonzero((local[:, 0] <= A.hi[i]) & (local[:, 1] >= A.lo[i]))
            if len(c):
                an = np.stack([directions[j] for j in c])
                best = np.argsort(-(an @ u))[:8]
                v, e = nnls(an[best].T, u, maxiter=500)
                if e**2 < 1 - 0.98**2:
                    j = c[best[np.argmax(v)]]
                    if np.dot(an[best[np.argmax(v)]], u) > 0.98:
                        bands[j] = [
                            min(bands[j][0], A.lo[i]),
                            max(bands[j][1], A.hi[i]),
                        ]
                    continue
        # A detector's mean spectrum outside the reference window is unused.
        bands[len(P)] = A.lo[i], A.hi[i]
        P.append(p)
        directions.append(u)
    P = np.asarray(P).reshape(-1, q)
    bands = bands[:len(P)].copy()
    if not len(P):
        return P, [np.zeros((len(s.t), 0)) for s in data], bands
    return P, unmix_profiles(data, A, P, bands), bands


def unmix_profiles(data, A, P, bands, extracted=None):
    unit = P / np.maximum(np.linalg.norm(P, axis=1)[:, None], 1e-60)
    H = []
    for si, s in enumerate(data):
        if extracted is None:
            coef = np.zeros((len(s.S), len(P)))
            for ii in np.flatnonzero(A.sample == si):
                j = A.index[ii]
                c = np.flatnonzero(
                    (bands[:, 0] <= A.hi[ii]) & (bands[:, 1] >= A.lo[ii])
                )
                if not len(c):
                    continue
                p = A.p[ii]
                take = c[np.argsort(-(unit[c] @ p))[:8]]
                w, _ = nnls(
                    unit[take].T, p / max(np.linalg.norm(p), 1e-60), maxiter=500
                )
                allocation = w * unit[take].sum(1)
                if allocation.sum() > 0:
                    coef[j, take] = s.norms[j] * allocation / allocation.sum()
            h = np.asarray(s.h @ coef)
        else:
            h = extracted[si]
        # Use the PSF allocation to locate each candidate in this sample.
        mask = np.zeros_like(h, dtype=np.bool_)
        rr = (
            np.median(s.warp[s.grid[:, 1]] - s.warp[s.grid[:, 0]])
            if len(s.grid)
            else 0.0
        )
        for j in range(len(P)):
            v = h[:, j]
            i = np.flatnonzero(v >= 0.01 * v.max()) if v.max() > 0 else []
            if len(i):
                left = s.warp[i[0]]
                right = s.warp[i[-1]]
            elif extracted is not None:
                left, right = bands[j]
            else:
                continue
            a, b = np.searchsorted(s.warp, [left - rr, right + rr])
            mask[a:b, j] = True
        h[~mask] = 0.0
        h = scan_unmix(s.y, s.bg, P, np.ascontiguousarray(h), mask, s.phi)
        for j in np.flatnonzero(h.sum(0) > 0):
            old = h[:, j]
            nz = np.flatnonzero(old > 0)
            sl = slice(nz[0], nz[-1] + 1)
            v = log_shape(s.t[sl], old[sl])
            old[sl] = v * (old.sum() / max(v.sum(), 1e-30))
        H.append(h)
    return H
