from types import SimpleNamespace
import numpy as np
from scipy import sparse
from .blocks import profile_sample, update_spectra
from .kernels import face, log_shape, nnls_gram, solve_p
from .workers import ordered_map
from .objective import means, objective_terms, score
from .support import contract, support_block
from .exchange import exchange


def unit(x):
    return x / np.maximum(np.linalg.norm(x, axis=-1, keepdims=True), 1e-30)


def atoms(data):
    sample = []
    index = []
    lo = []
    hi = []
    peak = []
    spectra = []
    height = []
    for i, s in enumerate(data):
        s.norms = s.S.sum(1)
        a, b, c = s.grid.T
        spectra.append(s.S / s.norms[:, None])
        n = len(a)
        sample.extend([i] * n)
        index.extend(range(n))
        lo.extend(s.warp[a])
        hi.extend(s.warp[b])
        peak.extend(s.warp[c])
        if n:
            height.extend(np.asarray(s.h.max(axis=0).toarray()).ravel() * s.norms)
    return SimpleNamespace(
        sample=np.array(sample),
        index=np.array(index),
        lo=np.array(lo),
        hi=np.array(hi),
        peak=np.array(peak),
        p=np.vstack(spectra),
        height=np.array(height),
    )


def group_dictionary(data, A):
    # Star clustering prevents spectral chaining. Overlap may extend along a
    # broad peak while every member must match the same spectral direction.
    norm = unit(A.p)
    live = np.ones(len(A.p), bool)
    clusters = []
    for j in np.argsort(-A.height):
        if not live[j]:
            continue
        ok = live & (norm @ norm[j] >= 0.97)
        lo, hi = A.lo[j], A.hi[j]
        while True:
            group = np.flatnonzero(ok & (A.hi >= lo) & (A.lo <= hi))
            left = A.lo[group].min()
            right = A.hi[group].max()
            if left == lo and right == hi:
                break
            lo = left
            hi = right
        live[group] = False
        clusters.append(group)
    P = []
    H = [np.zeros((len(s.t), len(clusters))) for s in data]
    bands = []
    for k, grp in enumerate(clusters):
        # Spectrum sharing and overlapping PSF fusion are the initialisation.
        w = np.sqrt(A.height[grp])
        P.append(np.average(A.p[grp], axis=0, weights=w))
        bands.append([A.lo[grp].min(), A.hi[grp].max()])
        for si in np.unique(A.sample[grp]):
            s = data[si]
            idx = A.index[grp[A.sample[grp] == si]]
            v = np.asarray(s.h[:, idx] @ s.norms[idx]).ravel()
            nn = np.flatnonzero(v > 0)
            sl = slice(nn[0], nn[-1] + 1)
            h = log_shape(s.t[sl], v[sl])
            H[si][sl, k] = h * (v.sum() / max(h.sum(), 1e-30))
    return np.array(P).reshape(-1, A.p.shape[1]), H, np.array(bands).reshape(-1, 2)


def update_background(data, mu, track=False):
    def one(args):
        s, m = args
        gain = 0.0
        for j, h in enumerate(s.Z.T):
            ix = np.flatnonzero(h > 0)
            hs = h[ix]
            old = s.alpha[j].copy()
            o = np.maximum(1e-8, m[ix] - hs[:, None] * old)
            q, g = solve_p(np.asfortranarray(s.y[ix] / s.phi), np.asfortranarray(o / s.phi), hs / s.phi, old, 0.0, track)
            gain += g
            delta = hs[:, None] * (q - old)
            m[ix] += delta
            s.bg[ix] += delta
            s.alpha[j] = q
        np.maximum(s.bg, 1e-8, out=s.bg)
        np.maximum(m, 1e-8, out=m)
        return gain

    return sum(ordered_map(one, zip(data, mu)))


def spectral_block(data, H, P, mu, price, track=False):
    if len(data) > 1 and len(P):
        return update_spectra(data, H, P, mu, price, track) + update_background(data, mu, track)
    gain = 0.0
    mass = np.sum([h.sum(0) for h in H], axis=0)
    for k in np.argsort(-mass):
        if mass[k] <= 0:
            continue
        old = P[k].copy()
        blocks = []
        yy = []
        oo = []
        hh = []
        for si, (s, h, m) in enumerate(zip(data, H, mu)):
            ix = np.flatnonzero(h[:, k] > 0)
            if not len(ix):
                continue
            a = h[ix, k]
            o = np.maximum(1e-8, m[ix] - a[:, None] * old)
            yy.append(s.y[ix] / s.phi)
            oo.append(o / s.phi)
            hh.append(a / s.phi)
            blocks.append((si, ix, a, o))
        q, g = solve_p(
            np.vstack(yy), np.vstack(oo), np.concatenate(hh), old, price, track
        )
        gain += g
        total = q.sum()
        for si, ix, a, o in blocks:
            mu[si][ix] = o + a[:, None] * q
            H[si][ix, k] = a * total
        P[k] = q / total if total > 0 else 0.0
    return gain + update_background(data, mu, track)


def profile_block(data, H, P, mu, bands, price, allow_birth):
    order = np.argsort(-np.sum([h.sum(0) for h in H], axis=0))

    def one(args):
        s, h, m = args
        psf_width = np.median(s.warp[s.grid[:,1]] - s.warp[s.grid[:,0]]) if len(s.grid) else 0.
        width = max(psf_width, 3*np.median(np.diff(s.warp)))
        return profile_sample(s.t, s.warp, s.y, s.phi, h, m, P, order, bands, width, price, allow_birth)

    gains = ordered_map(one, zip(data, H, mu))
    return np.array(gains).T[order].sum()


def lift(data, H, P, g):
    """Original-channel WLS initialisation conditional on the coarse profiles.

    The following likelihood blocks refine this starting value. Per-sample
    spline terms are included in the same stacked design.
    """
    profiles = sparse.vstack([sparse.csr_matrix(h) for h in H], format="csr")
    bases = sparse.block_diag([sparse.csr_matrix(s.Z) for s in data], format="csr")
    X = sparse.hstack([profiles, bases], format="csc")
    Y = np.vstack([s.Y for s in data])
    weights = np.concatenate([1 / (s.phi * np.maximum(s.Y.sum(1), 1)) for s in data])
    W = X.multiply(np.sqrt(weights)[:, None]).tocsc()
    scale = np.sqrt(np.asarray(W.power(2).sum(0)).ravel())
    scale[scale == 0] = 1
    W = W @ sparse.diags(1 / scale)
    G = np.ascontiguousarray((W.T @ W).toarray())
    U = np.ascontiguousarray(W.T @ (Y * np.sqrt(weights)[:, None]))
    # Disaggregate each provisional spectrum only for the starting point.
    ion = Y.sum(0)
    fra = ion / np.maximum(np.bincount(g, weights=ion)[g], 1e-100)
    start = np.vstack([P[:, g] * fra, *[s.alpha[:, g] * fra for s in data]])
    B = nnls_gram(G, U, np.ascontiguousarray(start * scale[:, None])) / scale[:, None]
    K = len(P)
    Q = B[:K]
    sizes = Q.sum(1)
    Q = Q / np.maximum(sizes[:, None], 1e-100)
    idx = K
    for i, s in enumerate(data):
        H[i] *= sizes
        s.alpha = B[idx : idx + s.Z.shape[1]].copy()
        idx += s.Z.shape[1]
        s.y = s.Y
        s.bg = np.maximum(1e-8, s.Z @ s.alpha + 1e-8)
    return Q


def snapshot(record, label, cycle, block, data, H, P, mu, price, D):
    if record is None:
        return
    actual = score(data, H, P, mu, price)
    if not np.isclose(actual["D"], D, rtol=3e-8, atol=1e-5):
        raise ArithmeticError(
            ("incremental deviance", label, cycle, block, D, actual["D"])
        )
    record.append(
        dict(stage=label, cycle=cycle, block=block, channels=P.shape[1], **actual)
    )


def optimize(
    data,
    H,
    P,
    bands,
    price,
    limit,
    record,
    label,
    birth=False,
    fixed=False,
    active=False,
    carry=None,
):
    if carry is None:
        mu = means(data, H, P)
        previous = score(data, H, P, mu, price)
    else:
        mu, previous = carry
    D = previous["D"]
    snapshot(record, label, 0, "start", data, H, P, mu, price, D)
    for cycle in range(1, limit + 1):
        if not fixed:
            D -= profile_block(
                data, H, P, mu, bands, price, allow_birth=birth and cycle == 1
            )
            snapshot(record, label, cycle, "profiles", data, H, P, mu, price, D)
        D -= spectral_block(data, H, P, mu, price, track=True)
        snapshot(
            record, label, cycle, "spectra and background", data, H, P, mu, price, D
        )
        if active:
            D -= support_block(data, H, P, mu, price)
            if not fixed and cycle in (2, 5):
                D -= contract(data, H, P, mu, bands, price)
            snapshot(record, label, cycle, "active support", data, H, P, mu, price, D)
        current = objective_terms(data, H, P, D, price)
        change = previous["objective"] - current["objective"]
        if change < -1e-5 * max(1.0, previous["objective"]):
            raise ArithmeticError((label, cycle, previous, current))
        previous = current
        if change <= 1e-3 * (current["objective"] + change):
            break
    return mu, previous


def refine_support(data, H, P, bands, price, passes, record, carry):
    mu, current = carry
    for attempt in range(1, passes + 1):
        moves, gain = exchange(data, H, P, mu, bands, price, with_gain=True)
        updated = objective_terms(data, H, P, current["D"] - gain, price)
        if updated["objective"] > current["objective"] + 1e-4:
            raise ArithmeticError(("support exchange increased GIC", current, updated))
        snapshot(
            record,
            "support exchange",
            attempt,
            "exchange",
            data,
            H,
            P,
            mu,
            price,
            updated["D"],
        )
        current = updated
        if not moves:
            break
        mu, current = optimize(
            data,
            H,
            P,
            bands,
            price,
            10,
            record,
            f"conditional refit {attempt}",
            carry=(mu, current),
        )
    return mu, current


def trim(H, P, bands):
    keep = (np.sum([h.sum(0) for h in H], axis=0) > 0) & (P.sum(1) > 0)
    if keep.all():
        return H, P, bands
    return [h[:, keep] for h in H], P[keep], bands[keep]
