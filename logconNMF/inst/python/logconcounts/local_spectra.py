"""Joint profile and spectrum exchange on a local nonnegative spectral cone."""

import numpy as np
from scipy.optimize import nnls
from .kernels import face, log_shape, update_h, update_p
from .objective import means, cell_deviance_difference
from .exchange import choices


def propose(H, P, bands):
    # Prefer templates observed independently of the possible mixture.
    seen = np.array([h.sum(0) > 0 for h in H])
    n = seen.sum(0)
    unit = P / np.maximum(np.linalg.norm(P, axis=1)[:, None], 1e-80)
    offers = []
    for k in np.flatnonzero(n):
        overlap = (bands[:, 0] <= bands[k, 1]) & (bands[:, 1] >= bands[k, 0])
        c = np.flatnonzero(overlap & (n > n[k]) & seen[~seen[:, k]].any(0))
        for dest in choices(k, c, P, unit):
            ids = np.array(dest)
            a, err = nnls(P[ids].T, P[k])
            rel = err / np.linalg.norm(P[k])
            if np.all(a > 0):
                offers.append((rel, k, ids, a))
    return sorted(offers, key=lambda x: x[0])


def get_blocks(data, H, P, mu, ids):
    blocks = []
    for si, (s, h, m) in enumerate(zip(data, H, mu)):
        pos = np.flatnonzero(h[:, ids].sum(1) > 0)
        if not len(pos):
            continue
        t = np.arange(pos[0], pos[-1] + 1)
        v = h[np.ix_(t, ids)]
        old = m[t]
        off = np.maximum(s.bg[t], old - v @ P[ids])
        blocks.append((si, t, s.t[t], s.y[t], off, old, v, s.phi))
    return blocks


def rent(t, v):
    return sum(face(t, col) for col in v.T)


def spectrum_update(blocks, vs, ps, ms, price):
    for j in np.argsort(-sum(v.sum(0) for v in vs)):
        parts = []
        yy = []
        oo = []
        xx = []
        old = ps[j]
        for b, v, m in zip(blocks, vs, ms):
            si, rows, t, y, off, orig, base, phi = b
            at = np.flatnonzero(v[:, j] > 0)
            if not len(at):
                continue
            a = v[at, j]
            o = np.maximum(off[at], m[at] - a[:, None] * old)
            parts.append((v, m, at, a, o))
            yy.append(y[at] / phi)
            oo.append(o / phi)
            xx.append(a / phi)
        if not parts:
            ps[j] = 0
            continue
        q = update_p(np.vstack(yy), np.vstack(oo), np.concatenate(xx), old, price)
        size = q.sum()
        for v, m, at, a, o in parts:
            v[at, j] = a * size
            m[at] = o + a[:, None] * q
        ps[j] = q / size if size else 0.0


def trial(bs, P, weights, price, retain):
    ps = P.copy()
    r = np.maximum(ps[0] - weights @ ps[1:], 0.0) if retain else np.zeros_like(ps[0])
    mass = r.sum()
    ps[0] = r / mass if mass else 0.0
    vs = []
    ms = []
    oldrent = (P > 0).sum() - len(P)
    for si, rows, t, y, off, old, v0, phi in bs:
        v = v0.copy()
        for j, w in enumerate(weights, 1):
            if w <= 0 or not v0[:, 0].any():
                continue
            u = v0[:, j] + w * v0[:, 0]
            u_shape = log_shape(t, u)
            v[:, j] = u_shape * (u.sum() / max(u_shape.sum(), 1e-80))
        v[:, 0] = v0[:, 0] * mass
        vs.append(v)
        ms.append(off + v @ ps)
        oldrent += rent(t, v0)
    # Each trial must improve the complete native criterion. Keep the best iterate.
    best = None
    for step in range(2):
        for b, v, m in zip(bs, vs, ms):
            si, rows, t, y, off, old, v0, phi = b
            for j in np.argsort(-v.sum(0)):
                if not v[:, j].any():
                    continue
                ion = np.flatnonzero(ps[j] > 0)
                p = ps[j, ion]
                o = np.maximum(off[:, ion], m[:, ion] - v[:, j, None] * p)
                u, _ = update_h(t, y[:, ion], o, p, v[:, j], phi, price)
                v[:, j] = u
                m[:, ion] = o + u[:, None] * p
        spectrum_update(bs, vs, ps, ms, price)
        d = sum(
            cell_deviance_difference(b[3], m, b[5], b[7]).sum() for b, m in zip(bs, ms)
        )
        active = sum(v.sum(0) for v in vs) > 0
        df = (
            (ps[active] > 0).sum()
            - active.sum()
            + sum(rent(b[2], v) for b, v in zip(bs, vs))
        )
        val = d + price * (df - oldrent)
        if val < -1e-5 and (best is None or val < best[0]):
            best = (
                val,
                d,
                ps.copy(),
                [v.copy() for v in vs],
                [m.copy() for m in ms],
                step + 1,
            )
    return best


def local_spectra(data, H, P, bands, price):
    mu = means(data, H, P)
    events = []
    handled = set()
    for rel, k, rec, w in propose(H, P, bands):
        if k in handled or not P[k].any() or not np.all(P[rec].sum(1) > 0):
            continue
        ids = np.r_[k, rec]
        bs = get_blocks(data, H, P, mu, ids)
        w, err = nnls(P[rec].T, P[k])
        rel = err / np.linalg.norm(P[k])
        if rel > 0.32 or np.any(w <= 0):
            continue
        proposals = [
            (b, keep)
            for keep in (True, False)
            if (b := trial(bs, P[ids], w, price, keep)) is not None
        ]
        if not proposals:
            continue
        (cost, d, ps, vs, ms, step), keep = min(proposals, key=lambda pair: pair[0][0])
        P[ids] = ps
        handled.add(k)
        for b, v, m in zip(bs, vs, ms):
            si, t = b[:2]
            H[si][np.ix_(t, ids)] = v
            mu[si][t] = m
        events.append(
            dict(
                source=k + 1,
                receivers=" ".join(str(r + 1) for r in rec),
                residual=keep,
                relative_cone_error=rel,
                iterations=step,
                D=d,
                GIC=cost,
            )
        )
    return events
