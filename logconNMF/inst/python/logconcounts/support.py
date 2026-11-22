"""Local changes of rank under the same deviance and face-dimension criterion."""

import numpy as np
from .kernels import face, log_shape, profile_gain, update_h, update_p


def support_block(data, H, P, mu, price):
    gain = 0.0
    for k in np.argsort(np.sum([h.sum(0) for h in H], axis=0)):
        ix = np.flatnonzero(P[k] > 0)
        if not len(ix):
            continue
        p = P[k, ix]
        merits = []
        positions = []
        deviances = []
        for i, (s, h, m) in enumerate(zip(data, H, mu)):
            t = np.flatnonzero(h[:, k] > 0)
            if not len(t):
                continue
            a = h[t, k]
            sl = np.ix_(t, ix)
            o = np.maximum(m[sl] - a[:, None] * p, 1e-8)
            improvement = profile_gain(s.y[sl], o, p, a, s.phi)
            deviances.append(improvement)
            merits.append(improvement - price * face(s.t[t], a))
            positions.append((i, t))
        empty = sum(max(0.0, x) for x in merits) <= price * (len(ix) - 1)
        for val, improvement, (i, t) in zip(merits, deviances, positions):
            if empty or val < 0:
                sl = np.ix_(t, ix)
                mu[i][sl] = np.maximum(mu[i][sl] - H[i][t, k, None] * p, 1e-8)
                H[i][t, k] = 0.0
                gain -= improvement
        if empty:
            P[k, :] = 0.0
    return gain


def contract_pair(data, H, P, mu, k, j, price):
    """Conditional comparison of the current two factors and one shared factor."""
    ions = np.flatnonzero((P[k] + P[j]) > 0)
    if not len(ions):
        return None
    a, b = P[k, ions], P[j, ions]
    blocks = []
    weights = np.zeros(2)
    old = -price * ((a > 0).sum() + (b > 0).sum() - 2)
    old_gain = 0.0
    for si, (s, h, m) in enumerate(zip(data, H, mu)):
        raw = h[:, k] + h[:, j]
        non = np.flatnonzero(raw > 0)
        if not len(non):
            continue
        t = np.arange(non[0], non[-1] + 1)
        sl = np.ix_(t, ions)
        ar = h[t, k]
        br = h[t, j]
        term = ar[:, None] * a + br[:, None] * b
        o = np.maximum(1e-8, m[sl] - term)
        y = s.y[sl]
        g = 2 / s.phi * (np.sum(y * np.log1p(term / o)) - term.sum())
        old_gain += g
        old += g - price * (face(s.t[t], ar) + face(s.t[t], br))
        v = log_shape(s.t[t], raw[t])
        v *= raw[t].sum() / max(1e-30, v.sum())
        blocks.append([si, t, y, o, v])
        weights += np.array([ar.sum(), br.sum()]) / s.phi
    if not blocks:
        return None
    p = (weights[0] * a + weights[1] * b) / weights.sum()
    for epoch in range(2):
        yy = []
        oo = []
        hh = []
        for bl in blocks:
            si, t, y, o, v = bl
            s = data[si]
            bl[4], _ = update_h(s.t[t], y, o, p, v, s.phi, price)
            live = bl[4] > 0
            if np.any(live):
                yy.append(y[live] / s.phi)
                oo.append(o[live] / s.phi)
                hh.append(bl[4][live] / s.phi)
        if not hh:
            return None
        q = update_p(np.vstack(yy), np.vstack(oo), np.concatenate(hh), p, price)
        if q.sum() <= 0:
            return None
        for bl in blocks:
            bl[4] *= q.sum()
        p = q / q.sum()
    new = -price * ((p > 0).sum() - 1)
    new_gain = 0.0
    for si, t, y, o, v in blocks:
        s = data[si]
        g = profile_gain(y, o, p, v, s.phi)
        new_gain += g
        new += g - price * face(s.t[t], v)
    if new <= old + 1e-6:
        return None
    P[j, :] = 0.0
    P[k, :] = 0.0
    P[k, ions] = p
    for si, t, y, o, v in blocks:
        H[si][t, k] = v
        H[si][t, j] = 0.0
        mu[si][np.ix_(t, ions)] = o + v[:, None] * p
    return new_gain - old_gain


def contract(data, H, P, mu, bands, price, max_neighbors=2):
    alive = P.sum(1) > 0
    mass = np.sum([h.sum(0) for h in H], axis=0)
    margin = np.median(
        [
            np.median(s.warp[s.grid[:, 1]] - s.warp[s.grid[:, 0]])
            for s in data
            if len(s.grid)
        ]
    )
    gain = 0.0
    order = np.argsort(-mass)
    for k in order:
        if not alive[k]:
            continue
        cand = np.flatnonzero(
            alive
            & (bands[:, 0] <= bands[k, 1] + margin)
            & (bands[:, 1] >= bands[k, 0] - margin)
        )
        cand = cand[(cand != k) & (mass[cand] <= mass[k])]
        if not len(cand):
            continue
        cos = (P[cand] @ P[k]) / np.maximum(
            np.linalg.norm(P[cand], axis=1) * np.linalg.norm(P[k]), 1e-60
        )
        for j in cand[np.argsort(-cos)[:max_neighbors]]:
            delta = contract_pair(data, H, P, mu, k, j, price)
            if delta is not None:
                alive[j] = False
                bands[k] = [
                    min(bands[k, 0], bands[j, 0]),
                    max(bands[k, 1], bands[j, 1]),
                ]
                mass[k] += mass[j]
                gain += delta
    return gain
