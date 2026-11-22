import numpy as np
from numba import njit


@njit(cache=True)
def knot_fit(t, y, w, knots):
    n = len(t)
    p = len(knots)
    diag = np.zeros(p)
    upper = np.zeros(p)
    rhs = np.zeros(p)
    k = 0
    for i in range(n):
        if k + 2 < p and i >= knots[k + 1]:
            k += 1
        b = (t[i] - t[knots[k]]) / (t[knots[k + 1]] - t[knots[k]])
        a = 1 - b
        diag[k] += w[i] * a * a
        diag[k + 1] += w[i] * b * b
        upper[k] += w[i] * a * b
        rhs[k] += w[i] * a * y[i]
        rhs[k + 1] += w[i] * b * y[i]
    for k in range(1, p):
        ratio = upper[k - 1] / diag[k - 1]
        diag[k] -= ratio * upper[k - 1]
        rhs[k] -= ratio * rhs[k - 1]
    rhs[-1] /= diag[-1]
    for k in range(p - 2, -1, -1):
        rhs[k] = (rhs[k] - upper[k] * rhs[k + 1]) / diag[k]
    return np.interp(t, t[knots], rhs)


@njit(cache=True)
def project(t, y, w):
    n = len(t)
    if n <= 2:
        return y.copy()
    # Remove the affine null space implicitly by fitting the two endpoints.
    x = (t - t[0]) / (t[-1] - t[0])
    w = w / w.max()
    active = np.zeros(n, np.bool_)
    active[0] = True
    active[-1] = True
    z = knot_fit(x, y, w, np.flatnonzero(active))
    for it in range(10 * n + 10):
        knots = np.flatnonzero(active)
        v = knot_fit(x, y, w, knots)
        alpha = 1.0
        drop = -1
        for j in range(1, len(knots) - 1):
            a, b, c = knots[j - 1], knots[j], knots[j + 1]
            frac = (x[b] - x[a]) / (x[c] - x[a])
            slack = v[b] - ((1 - frac) * v[a] + frac * v[c])
            if slack < -1e-10:
                old = max(0.0, z[b] - ((1 - frac) * z[a] + frac * z[c]))
                step = old / (old - slack)
                if step < alpha:
                    alpha = step
                    drop = b
        if drop >= 0:
            z += alpha * (v - z)
            active[drop] = False
            continue
        z = v
        # Direction -max(x-x_j,0) introduces one admissible downward bend.
        g = w * (z - y)
        total = 0.0
        moment = 0.0
        best = 1e-12 * max(1.0, np.abs(g).sum())
        add = -1
        for i in range(n - 1, 0, -1):
            total += g[i]
            moment += x[i] * g[i]
            k = i - 1
            score = moment - x[k] * total
            if not active[k] and score > best:
                best = score
                add = k
        if add < 0:
            return z
        active[add] = True
    return z


@njit(cache=True)
def log_shape(t, h):
    """Project log ordinates on the concave cone; keep a contiguous support."""
    h = np.maximum(h, 0.0)
    ix = np.flatnonzero(h > 0)
    if not len(ix):
        return h
    a, b = ix[0], ix[-1] + 1
    v = h[a:b]
    hi = v.max()
    w = np.maximum(v / hi, 1e-6)
    z = project(t[a:b], np.log(np.maximum(v, hi * 1e-9)), w)
    out = np.zeros(len(h))
    out[a:b] = np.exp(np.minimum(z, np.log(hi) + 4))
    # Capping the logarithm at a constant preserves concavity.
    out[a:b] = np.maximum(out[a:b], 0.0)
    return out


@njit(cache=True)
def face(t, h):
    ix = np.flatnonzero(h > 0)
    if len(ix) < 3:
        return len(ix)
    z = np.log(h[ix])
    tol = 1.4901161193847656e-8 * (1 + np.max(np.abs(z)))
    d = 2
    for j in range(1, len(ix) - 1):
        a, b, c = ix[j - 1], ix[j], ix[j + 1]
        gap = z[j] - z[j - 1] - (z[j + 1] - z[j - 1]) * (t[b] - t[a]) / (t[c] - t[a])
        if gap > tol:
            d += 1
    return d


@njit(cache=True)
def curve_error(t, h):
    i = np.flatnonzero(h > 0)
    if len(i) and i[-1] - i[0] + 1 != len(i):
        return np.inf
    if len(i) < 3:
        return 0.0
    z = np.log(h[i])
    t = t[i]
    return max(
        0.0,
        np.max(
            z[:-2] + (z[2:] - z[:-2]) * (t[1:-1] - t[:-2]) / (t[2:] - t[:-2]) - z[1:-1]
        ),
    )


@njit(cache=True)
def profile_gain(y, o, p, h, phi):
    total = 0.0
    for i in range(len(h)):
        if h[i] <= 0:
            continue
        for j in range(len(p)):
            v = h[i] * p[j]
            if v > 0:
                total += y[i, j] * np.log1p(v / o[i, j]) - v
    return 2 * total / phi


@njit(cache=True)
def rescale(y, o, p, h):
    hp = h.max()
    if hp <= 0:
        return h.copy()
    # Safeguarded Newton in peak-height units.
    hh = h / hp
    norm = hh.sum() * p.sum()
    lo = 0.0
    hi = y.sum() / max(norm, 1e-30)
    a = min(hi, hp)
    for _ in range(18):
        g = -norm
        c = 0.0
        for i in range(len(h)):
            if hh[i] == 0:
                continue
            for j in range(len(p)):
                v = hh[i] * p[j]
                if v == 0:
                    continue
                z = v / (o[i, j] + a * v)
                g += y[i, j] * z
                c += y[i, j] * z * z
        if abs(g) <= 1e-8 * norm:
            break
        if g > 0:
            lo = a
        else:
            hi = a
        if hi <= hp * 1e-12:
            return np.zeros(len(h))
        new = a + g / max(c, 1e-60)
        a = new if lo < new < hi else 0.5 * (lo + hi)
    return hh * a


@njit(cache=True)
def scan_coefficients(y, o, p, h):
    n, m = y.shape
    v = np.zeros(n)
    gain = np.zeros(n)
    info = np.zeros(n)
    sump = p.sum()
    for i in range(n):
        g0 = -sump
        for j in range(m):
            g0 += p[j] * y[i, j] / o[i, j]
        if g0 <= 0:
            continue
        lo = 0.0
        hi = y[i].sum() / sump
        a = min(hi, max(h[i], 0.001 * hi))
        for _ in range(19):
            g = -sump
            c = 0.0
            for j in range(m):
                z = p[j] / (o[i, j] + a * p[j])
                g += y[i, j] * z
                c += y[i, j] * z * z
            if abs(g) < 1e-7 * sump:
                break
            if g > 0:
                lo = a
            else:
                hi = a
            new = a + g / max(c, 1e-60)
            a = new if lo < new < hi else 0.5 * (lo + hi)
        v[i] = a
        for j in range(m):
            term = a * p[j]
            gain[i] += y[i, j] * np.log1p(term / o[i, j]) - term
            info[i] += p[j] * p[j] / (o[i, j] + term)
    return v, gain, info


@njit(cache=True)
def update_h(t, y, o, p, h, phi, price):
    old = profile_gain(y, o, p, h, phi) - price * face(t, h)
    best = h.copy()
    best_score = old
    if best_score < 0:
        best[:] = 0
        best_score = 0.0
    v, g, info = scan_coefficients(y, o, p, h)
    # Best contiguous interval for an initially saturated scan-rate proposal.
    # The actual comparison below pays only for free log-shape ordinates.
    right = 0
    left = 0
    start = 0
    cur = 0.0
    value = 0.0
    for i in range(len(g)):
        cur += 2 * g[i] / phi - price
        if cur > value:
            value = cur
            left = start
            right = i + 1
        if cur <= 0:
            start = i + 1
            cur = 0.0
    if right > left:
        left = max(0, left - 1)
        right = min(len(t), right + 1)
        v0 = v[left:right]
        mx = v0.max()
        if mx > 0:
            w = np.maximum(info[left:right] * v0 * v0, 1e-15)
            w = np.maximum(w, w.max() * 1e-7)
            z = project(t[left:right], np.log(np.maximum(v0, mx * 1e-8)), w)
            new = np.zeros(len(h))
            new[left:right] = np.exp(np.minimum(z, np.log(mx) + 3))
            new = rescale(y, o, p, new)
            val = profile_gain(y, o, p, new, phi) - price * face(t, new)
            if val > best_score + 1e-7:
                best_score = val
                best = new
    if h.max() > 0:
        scaled = rescale(y, o, p, h)
        val = profile_gain(y, o, p, scaled, phi) - price * face(t, scaled)
        if val > best_score + 1e-7:
            best_score = val
            best = scaled
    return best, best_score - old


@njit(cache=True)
def solve_p(y, o, h, old, price, track_gain):
    # The inputs y, o and h have already been divided by sample dispersion.
    n, m = y.shape
    v = np.zeros(m)
    total = h.sum()
    improvement = 0.0
    if total <= 0:
        return v, improvement
    for j in range(m):
        if track_gain and old[j] > 0:
            improvement += old[j] * total
            for i in range(n):
                improvement -= y[i, j] * np.log1p(old[j] * h[i] / o[i, j])
        der = -total
        for i in range(n):
            der += y[i, j] * h[i] / o[i, j]
        if der <= 0:
            continue
        lo = 0.0
        hi = y[:, j].sum() / total
        a = min(hi, max(old[j], 1e-4 * hi))
        for _ in range(25):
            der = -total
            c = 0.0
            for i in range(n):
                r = h[i] / (o[i, j] + a * h[i])
                der += y[i, j] * r
                c += y[i, j] * r * r
            if abs(der) < 1e-8 * total:
                break
            if der > 0:
                lo = a
            else:
                hi = a
            new = a + der / max(c, 1e-90)
            a = new if lo < new < hi else 0.5 * (lo + hi)
        gain = -a * total
        for i in range(n):
            gain += y[i, j] * np.log1p(a * h[i] / o[i, j])
        if 2 * gain > price:
            v[j] = a
            if track_gain:
                improvement += gain
    return v, 2 * improvement


@njit(cache=True)
def update_p(y, o, h, old, price):
    return solve_p(y, o, h, old, price, False)[0]


@njit(cache=True)
def nnls_gram(G, U, B):
    """Many nonnegative right-hand sides with one positive-semidefinite Gram."""
    B = B.copy()
    R = U - G @ B
    for _ in range(120):
        mx = 0.0
        for j in range(len(B)):
            if G[j, j] <= 0:
                continue
            change = np.maximum(0.0, B[j] + R[j] / G[j, j]) - B[j]
            mx = max(mx, np.max(np.abs(change)))
            B[j] += change
            for k in np.flatnonzero(G[:, j] != 0):
                R[k] -= G[k, j] * change
        if mx < 1e-8 * max(1.0, np.max(B)):
            break
    return B
