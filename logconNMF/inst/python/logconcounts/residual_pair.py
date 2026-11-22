"""Test two shared spectral directions when samples disagree about one direction."""

import numpy as np
from .anchors import scan_unmix
from .kernels import face, log_shape, rescale, update_h, update_p
from .objective import means, cell_deviance_difference


def occurrences(data, H, P, mu, k):
    out = []
    for si, (s, h, m) in enumerate(zip(data, H, mu)):
        rows = np.flatnonzero(h[:, k] > 0)
        if not len(rows):
            continue
        v = h[rows, k]
        o = np.maximum(s.bg[rows], m[rows] - v[:, None] * P[k])
        out.append((si, rows, v, s.y[rows], o, m[rows]))
    return out


def partition(data, obs, p):
    if len(obs) < 2:
        return None
    G = []
    I = []
    for si, rows, h, y, o, m in obs:
        phi = data[si].phi
        H = h[:, None]
        W = 1 / (m * phi)
        G.append(((y - m) * W * H).sum(0))
        I.append((H * H * W).sum(0))
    g = np.array(G)
    w = np.maximum(I, 1e-80)
    # Remove the scale direction; each fitted occurrence already has its own area.
    coef = (g * p).sum(1) / np.maximum((w * p * p).sum(1), 1e-80)
    g -= coef[:, None] * w * p
    f = np.maximum(0.0, p + g / w)
    f /= np.maximum(f.sum(1)[:, None], 1e-60)
    root = np.sqrt(f)
    _, _, vt = np.linalg.svd(root - root.mean(0), full_matrices=False)
    order = np.argsort(root @ vt[0])

    def gain(gg, ii):
        q = np.maximum(-p, gg / ii)
        return (2 * gg * q - ii * q * q).sum()

    agg = g.sum(0)
    info = w.sum(0)
    base = gain(agg, info)
    # Enumerate all cuts on the leading spectral contrast and all singletons.
    cuts = [order[:i] for i in range(1, len(obs))] + [
        np.array([i]) for i in range(len(obs))
    ]
    best = -np.inf
    sel = None
    for ids in cuts:
        take = np.zeros(len(obs), bool)
        take[ids] = True
        gg = g[take].sum(0)
        ii = w[take].sum(0)
        score = gain(gg, ii) + gain(g[~take].sum(0), w[~take].sum(0)) - base
        if score > best:
            sel = take
            best = score
    return best, sel


def split_offer(data, obs, p, mask, price):
    groups = []
    dd = 0.0
    df = -(np.count_nonzero(p) - 1)
    for ii in (np.flatnonzero(mask), np.flatnonzero(~mask)):
        rows = [obs[i] for i in ii]
        ys = np.vstack([y / data[s].phi for s, t, h, y, o, m in rows])
        os = np.vstack([o / data[s].phi for s, t, h, y, o, m in rows])
        hs = np.concatenate([h / data[s].phi for s, t, h, y, o, m in rows])
        q = update_p(ys, os, hs, p, price)
        scale = q.sum()
        if not scale:
            return None
        q /= scale
        df += np.count_nonzero(q) - 1
        values = []
        ions = np.flatnonzero(q > 0)
        for si, t, h, y, o, m in rows:
            s = data[si]
            v, _ = update_h(
                s.t[t], y[:, ions], o[:, ions], q[ions], h * scale, s.phi, price
            )
            new = o + v[:, None] * q
            dd += cell_deviance_difference(y, new, m, s.phi).sum()
            df += face(s.t[t], v) - face(s.t[t], h)
            values.append((si, t, v, new))
        if not any(v.any() for si, t, v, m in values):
            return None
        groups.append((q, values))
    return dd + price * df, dd, groups


def replicate_gain(data, obs, p, mask, price):
    """Cross-check spectral heterogeneity using one held-out occurrence at a time."""
    group = [np.flatnonzero(mask), np.flatnonzero(~mask)]
    if min(map(len, group)) < 2:
        return False

    def train(ids):
        rows = [obs[i] for i in ids]
        q = update_p(
            np.vstack([y / data[s].phi for s, t, h, y, o, m in rows]),
            np.vstack([o / data[s].phi for s, t, h, y, o, m in rows]),
            np.concatenate([h / data[s].phi for s, t, h, y, o, m in rows]),
            p,
            price,
        )
        return q / max(q.sum(), 1e-80)

    for members in group:
        score = 0.0
        for i in members:
            si, t, h, y, o, m = obs[i]
            s = data[si]
            local = train(members[members != i])
            pooled = train(np.delete(np.arange(len(obs)), i))
            if not local.any() or not pooled.any():
                return False
            # Profile out the nuisance area under each held-out prediction.
            v = rescale(y, o, local, h)
            u = rescale(y, o, pooled, h)
            score += cell_deviance_difference(
                y, o + u[:, None] * pooled, o + v[:, None] * local, s.phi
            ).sum()
        if score <= 0:
            return False
    return True


def admit_coelution(data, H, mu, ids, p, obs, price):
    gain = dev = 0.0
    moves = 0
    ion = np.flatnonzero(p.sum(0) > 0)
    p = p[:, ion]
    for si, t, h, y, off, old in obs:
        s = data[si]
        y = y[:, ion]
        o = off[:, ion]
        at = np.ix_(t, ion)
        old = mu[si][at]
        prev = H[si][np.ix_(t, ids)]
        v = scan_unmix(y, o, p, prev.copy(), np.ones_like(prev, dtype=bool), s.phi)
        for j in range(2):
            z = log_shape(s.t[t], v[:, j])
            v[:, j] = z * (v[:, j].sum() / max(z.sum(), 1e-80))
        m = o + v @ p
        for j in np.argsort(-v.sum(0)):
            use = np.flatnonzero(p[j] > 0)
            q = p[j, use]
            rest = np.maximum(o[:, use], m[:, use] - v[:, j, None] * q)
            v[:, j], _ = update_h(s.t[t], y[:, use], rest, q, v[:, j], s.phi, price)
            m[:, use] = rest + v[:, j, None] * q
        dd = cell_deviance_difference(y, m, old, s.phi).sum()
        dj = dd + price * sum(
            face(s.t[t], x) - face(s.t[t], z) for x, z in zip(v.T, prev.T)
        )
        if dj < -1e-5:
            mu[si][at] = m
            H[si][np.ix_(t, ids)] = v
            dev += dd
            gain += dj
            moves += 1
    return gain, dev, moves


def residual_pair(
    data, H, P, bands, price, replicated=False, coexisting=False, polished=False
):
    mu = means(data, H, P)
    op = []
    events = []
    for k, p in enumerate(P):
        if not p.any():
            continue
        obs = occurrences(data, H, P, mu, k)
        z = partition(data, obs, p)
        if z is None:
            continue
        gain, _ = z
        rent = price * (np.count_nonzero(p) - 1)
        if gain > rent:
            op.append((gain - rent, k))
    for _, k in sorted(op, reverse=True):
        p = P[k]
        obs = occurrences(data, H, P, mu, k)
        gain, mask = partition(data, obs, p)
        if gain <= price * (np.count_nonzero(p) - 1):
            continue
        if replicated and not replicate_gain(data, obs, p, mask, price):
            continue
        z = split_offer(data, obs, p, mask, price)
        if z is None or z[0] >= -1e-5:
            continue
        val, dd, groups = z
        j = len(P)
        P = np.vstack((P, groups[1][0]))
        P[k] = groups[0][0]
        bands = np.vstack((bands, bands[k]))
        for si, h in enumerate(H):
            H[si] = np.column_stack((h, np.zeros(len(h))))
        for new, pair in zip((k, j), groups):
            for si, t, v, m in pair[1]:
                H[si][t, k] = 0.0
                H[si][t, new] = v
                mu[si][t] = m
        extra = 0
        if coexisting:
            delta, ddev, extra = admit_coelution(
                data, H, mu, [k, j], P[[k, j]], obs, price
            )
            val += delta
            dd += ddev
        if polished and extra:
            from .pair_polish import polish

            delta, ddev = polish(data, H, P, mu, [k, j], price)
            val += delta
            dd += ddev
        events.append(
            dict(
                source=k + 1,
                new_id=j + 1,
                samples_a=len(groups[0][1]),
                samples_b=len(groups[1][1]),
                predicted_gain=gain,
                D=dd,
                GIC=val,
                joint_windows=extra,
            )
        )
    return P, bands, events
