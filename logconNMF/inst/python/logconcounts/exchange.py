"""One joint active-set move for a spectral ray represented by nearby rays."""

import numpy as np
from .kernels import face, log_shape, update_h
from .anchors import scan_unmix


def choices(k, c, P, norm):
    if not len(c):
        return []
    u = norm[k]
    x = norm[c]
    cos = x @ u
    top = np.argsort(-cos)[:8]
    out = []
    if cos[top[0]] > 0.90:
        out.append((c[top[0]],))
    best = 0.90**2
    pair = None
    for ai, i in enumerate(top):
        for j in top[ai + 1 :]:
            cross = x[i] @ x[j]
            det = 1 - cross * cross
            if det <= 1e-10:
                continue
            a = (cos[i] - cross * cos[j]) / det
            b = (cos[j] - cross * cos[i]) / det
            if min(a, b) <= 0:
                continue
            r = a * cos[i] + b * cos[j]
            if r > best:
                best = r
                pair = (c[i], c[j])
    if pair:
        out.append(pair)
    return out


def trial(data, H, P, mu, si, k, rec, price):
    s = data[si]
    h = H[si]
    ids = [k, *rec]
    keep = np.flatnonzero(P[ids].sum(0) > 0)
    use = np.flatnonzero(h[:, ids].sum(1) > 0)
    t = np.arange(use[0], use[-1] + 1)
    sel = np.ix_(t, keep)
    p = np.ascontiguousarray(P[rec, :][:, keep])
    prior = h[np.ix_(t, ids)]
    oldmu = mu[si][sel]
    o = np.maximum(oldmu - prior @ P[ids, :][:, keep], 1e-8)
    y = s.y[sel]
    v = np.ascontiguousarray(prior[:, 1:].copy())
    v += prior[:, :1] / len(rec)
    v = scan_unmix(y, o, p, v, np.ones_like(v, dtype=bool), s.phi)
    for j in range(len(rec)):
        a = v[:, j]
        b = log_shape(s.t[t], a)
        v[:, j] = b * (a.sum() / max(b.sum(), 1e-60))
    m = o + v @ p
    for j in np.argsort(-v.sum(0)):
        pc = p[j]
        ix = pc > 0
        if not np.any(ix):
            continue
        off = np.maximum(m[:, ix] - v[:, j, None] * pc[ix], 1e-8)
        a, _ = update_h(s.t[t], y[:, ix], off, pc[ix], v[:, j], s.phi, price)
        m[:, ix] = off + a[:, None] * pc[ix]
        v[:, j] = a
    diff = m - oldmu
    rel = diff / oldmu
    lr = np.log(m) - np.log(oldmu)
    small = np.abs(rel) < 0.5
    lr[small] = np.log1p(rel[small])
    dd = 2 / s.phi * np.sum(diff - y * lr)
    d = dd + price * (
        sum(face(s.t[t], v[:, j]) for j in range(len(rec)))
        - sum(face(s.t[t], prior[:, j]) for j in range(len(ids)))
    )
    return d, (si, t, keep, rec, v, m), dd


def exchange(data, H, P, mu, bands, price, with_gain=False):
    mass = np.stack([h.sum(0) for h in H])
    norm = P / np.maximum(np.linalg.norm(P, axis=1)[:, None], 1e-60)
    width = np.median(
        [
            np.median(s.warp[s.grid[:, 1]] - s.warp[s.grid[:, 0]])
            for s in data
            if len(s.grid)
        ]
    )
    order = np.lexsort((mass.sum(0), (mass > 0).sum(0)))
    count = 0
    gain = 0.0
    for k in order:
        members = np.flatnonzero(mass[:, k] > 0)
        if not len(members):
            continue
        c = np.flatnonzero(
            (mass.sum(0) > 0)
            & (bands[:, 0] <= bands[k, 1] + width)
            & (bands[:, 1] >= bands[k, 0] - width)
        )
        c = c[c != k]
        opts = choices(k, c, P, norm)
        if not opts:
            continue
        offers = []
        for si in members:
            offers.append(
                min(
                    (trial(data, H, P, mu, si, k, list(rr), price) for rr in opts),
                    key=lambda a: a[0],
                )
            )
        retire = sum(o[0] for o in offers) < price * ((P[k] > 0).sum() - 1) - 1e-5
        for delta, (si, t, keep, rec, v, m), dd in offers:
            if not retire and delta >= -1e-5:
                continue
            H[si][t, k] = 0.0
            mass[si, k] = 0.0
            mu[si][np.ix_(t, keep)] = m
            for j, r in enumerate(rec):
                H[si][t, r] = v[:, j]
                mass[si, r] = v[:, j].sum()
                bands[r] = [
                    min(bands[r, 0], bands[k, 0]),
                    max(bands[r, 1], bands[k, 1]),
                ]
            count += 1
            gain -= dd
        if not np.any(mass[:, k] > 0):
            P[k, :] = 0.0
    return (count, gain) if with_gain else count
