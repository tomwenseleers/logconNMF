"""Refine two nearly parallel spectra using singly occupied windows for initialization."""

import numpy as np
from .kernels import log_shape, update_h, update_p
from .local_spectra import get_blocks, rent, spectrum_update
from .objective import cell_deviance_difference


def two_ray(y, o, p, start):
    """Fisher scoring for a two-variable nonnegative fit, exactly solving each quadratic."""
    h = start.copy()
    for it in range(4):
        m = np.maximum(o + h @ p, 1e-8)
        w = 1 / m
        a = w @ (p[0] * p[0])
        b = w @ (p[0] * p[1])
        c = w @ (p[1] * p[1])
        u = ((y - o) * w) @ p.T
        det = a * c - b * b
        stable = det > 1e-12 * a * c
        v = np.zeros_like(h)
        v[stable, 0] = (c[stable] * u[stable, 0] - b[stable] * u[stable, 1]) / det[
            stable
        ]
        v[stable, 1] = (a[stable] * u[stable, 1] - b[stable] * u[stable, 0]) / det[
            stable
        ]
        edge = ~stable | (v.min(1) < 0)
        single = np.maximum(u, 0.0) / np.maximum(np.column_stack((a, c)), 1e-100)
        pick = np.argmax(single * u, axis=1)
        v[edge] = 0.0
        ii = np.flatnonzero(edge)
        v[ii, pick[ii]] = single[ii, pick[ii]]
        h = v
    return h


def polish(data, H, P, mu, ids, price):
    bs = get_blocks(data, H, P, mu, ids)
    seen = np.array([b[6].sum(0) > 0 for b in bs])
    if not seen.all(1).any() or any(
        not (seen[:, j] & ~seen[:, 1 - j]).any() for j in range(2)
    ):
        return 0.0, 0.0
    old = P[ids].copy()
    ps = old.copy()
    # These singly occupied occurrences propose spectra. All original occurrences
    # and all ions participate in the acceptance test and shared-spectrum update.
    for j in range(2):
        obs = [b for b, s in zip(bs, seen) if s[j] and not s[1 - j]]
        q = update_p(
            np.vstack([b[3] / b[7] for b in obs]),
            np.vstack([b[4] / b[7] for b in obs]),
            np.concatenate([b[6][:, j] / b[7] for b in obs]),
            ps[j],
            price,
        )
        if not q.sum():
            return 0.0, 0.0
        ps[j] = q / q.sum()
    vs = [b[6].copy() for b in bs]
    best = None
    original_rent = (old > 0).sum() - 2 + sum(rent(b[2], b[6]) for b in bs)
    for it in range(2):
        ms = []
        for i, (si, rows, t, y, o, prior, v0, phi) in enumerate(bs):
            v = vs[i]
            ids0 = np.flatnonzero(seen[i])
            if len(ids0) == 2:
                v = two_ray(y, o, ps, v)
                for j in range(2):
                    h = log_shape(t, v[:, j])
                    v[:, j] = h * (v[:, j].sum() / max(h.sum(), 1e-100))
            m = o + v @ ps
            for j in ids0:
                ix = np.flatnonzero(ps[j] > 0)
                p = ps[j, ix]
                off = np.maximum(o[:, ix], m[:, ix] - v[:, j, None] * p)
                v[:, j], _ = update_h(t, y[:, ix], off, p, v[:, j], phi, price)
                m[:, ix] = off + v[:, j, None] * p
            vs[i] = v
            ms.append(m)
        spectrum_update(bs, vs, ps, ms, price)
        if not (ps.sum(1) > 0).all():
            break
        d = sum(
            cell_deviance_difference(b[3], m, b[5], b[7]).sum() for b, m in zip(bs, ms)
        )
        df = (
            (ps > 0).sum()
            - 2
            + sum(rent(b[2], v) for b, v in zip(bs, vs))
            - original_rent
        )
        cost = d + price * df
        if cost < -1e-5 and (best is None or cost < best[0]):
            best = (cost, d, ps.copy(), [v.copy() for v in vs], [m.copy() for m in ms])
    if best is None:
        return 0.0, 0.0
    dj, d, ps, vs, ms = best
    P[ids] = ps
    for b, v, m in zip(bs, vs, ms):
        si, t = b[:2]
        H[si][np.ix_(t, ids)] = v
        mu[si][t] = m
    return dj, d
