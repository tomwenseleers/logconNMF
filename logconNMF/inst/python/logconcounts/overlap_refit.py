"""GIC-checked reallocation of mixed, coeluting factors in an original-ion count fit.

data contains the observed ion matrix y, fitted background bg, native time t,
aligned time warp and quasi-Poisson dispersion phi for each sample. H holds
TIC profiles (scans x shared identities), and rows of P are unit-sum spectra.
H and P are updated in place. Zero columns are compacted by the caller.
"""

import numpy as np
from scipy.optimize import nnls
from .kernels import face, log_shape, update_h
from .local_spectra import MAX_SPECTRAL_RESIDUAL, get_blocks, trial
from .objective import cell_deviance_difference, means


def _aligned_ranges(data, H, ncomp):
    bounds = np.column_stack((np.full(ncomp, np.inf), np.full(ncomp, -np.inf)))
    for s, h in zip(data, H):
        for k in np.flatnonzero(h.max(0) > 0):
            ix = np.flatnonzero(h[:, k] > 0)
            bounds[k] = min(bounds[k, 0], s.warp[ix[0]]), max(
                bounds[k, 1], s.warp[ix[-1]]
            )
    return bounds


def _shared_exchanges(data, H, P, mu, price):
    bounds = _aligned_ranges(data, H, len(P))
    mass = np.array([h.sum(0) for h in H])
    ns = (mass > 0).sum(0)

    def representation(k):
        seen = mass[:, k] > 0
        if not seen.any() or not P[k].any():
            return None
        candidates = np.flatnonzero(
            (bounds[:, 0] <= bounds[k, 1])
            & (bounds[:, 1] >= bounds[k, 0])
            & (ns > ns[k])
            & (mass[~seen] > 0).any(0)
        )
        if not candidates.size:
            return None
        try:
            w, err = nnls(P[candidates].T, P[k], maxiter=1000)
        except (RuntimeError, np.linalg.LinAlgError):
            return None
        return err / np.linalg.norm(P[k]), candidates[w > 0], w[w > 0]

    offers = []
    for k in range(len(P)):
        z = representation(k)
        if z is not None and z[0] <= MAX_SPECTRAL_RESIDUAL:
            offers.append((z[0], k))
    moves = []
    for _, k in sorted(offers):
        z = representation(k)
        if z is None or z[0] > MAX_SPECTRAL_RESIDUAL:
            continue
        rel, rec, w = z
        # Reduce negligible spectral directions before the local count fits.
        n_eff = sum(mass[i, k] / s.phi for i, s in enumerate(data))
        while len(w) > 1:
            keep = np.arange(len(w)) != np.argmin(w * np.linalg.norm(P[rec], axis=1))
            try:
                ww, err = nnls(P[rec[keep]].T, P[k])
            except (RuntimeError, np.linalg.LinAlgError):
                break
            new_rel = err / np.linalg.norm(P[k])
            if (new_rel**2 - rel**2) * max(1.0, n_eff) > price:
                break
            rec, w, rel = rec[keep][ww > 0], ww[ww > 0], new_rel
        if not len(w) or rel > MAX_SPECTRAL_RESIDUAL:
            continue
        ids = np.r_[k, rec]
        blocks = get_blocks(data, H, P, mu, ids)
        opts = [
            z
            for retain in (False, True)
            if (z := trial(blocks, P[ids], w, price, retain)) is not None
        ]
        if not opts:
            continue
        dj, dd, ps, vs, ms, _ = min(opts, key=lambda z: z[0])
        P[ids] = ps
        for block, v, m in zip(blocks, vs, ms):
            si, t = block[:2]
            H[si][np.ix_(t, ids)] = v
            mu[si][t] = m
            mass[si, ids] = v.sum(0)
        ns[:] = (mass > 0).sum(0)
        moves.append(
            dict(
                stage="shared mixture",
                source=k + 1,
                receivers=" ".join(str(j + 1) for j in rec),
                GIC=dj,
                D=dd,
            )
        )
    return moves


def _overlaps(h):
    pieces = []
    for k in np.flatnonzero(h.max(0) > 0):
        rows = np.flatnonzero(h[:, k] >= 0.5 * h[:, k].max())
        pieces.append((rows[0], rows[-1], k))
    groups = []
    for lo, hi, k in sorted(pieces):
        if not groups or lo > end:
            groups.append([k])
            end = hi
        else:
            groups[-1].append(k)
            end = max(end, hi)
    return [np.array(g) for g in groups if len(g) > 1]


def _refit_window(s, h, P, mu, ids, ns, price):
    block = h[:, ids]
    pos = np.flatnonzero(block.sum(1) > 0)
    rows = np.arange(pos[0], pos[-1] + 1)
    ions = np.flatnonzero(P[ids].sum(0) > 0)
    at = np.ix_(rows, ions)
    q = np.ascontiguousarray(P[np.ix_(ids, ions)])
    all_rows = np.arange(len(ids))
    candidates = []
    for j in all_rows:
        try:
            _, err = nnls(q[all_rows != j].T, q[j])
        except (RuntimeError, np.linalg.LinAlgError):
            continue
        if err <= MAX_SPECTRAL_RESIDUAL * np.linalg.norm(q[j]):
            candidates.append(j)
    if not candidates:
        return None
    t, y, old = s.t[rows], s.y[at], mu[at]
    h0 = block[rows]
    off = np.maximum(s.bg[at], old - h0 @ q)
    df0 = sum(face(t, v) for v in h0.T)
    last = (ns[ids] == 1) * (np.count_nonzero(q, axis=1) - 1)

    def refit(keep, prior):
        spectra = q[keep]
        v = np.ascontiguousarray(prior[:, keep])
        for i in range(len(t)):
            for _ in range(3):
                w = 1 / np.sqrt(np.maximum(off[i] + v[i] @ spectra, 1.0))
                try:
                    v[i] = nnls(
                        np.ascontiguousarray(w[:, None] * spectra.T),
                        w * (y[i] - off[i]),
                        maxiter=3000,
                    )[0]
                except (RuntimeError, np.linalg.LinAlgError):
                    return None
        for j in range(len(keep)):
            shape = log_shape(t, v[:, j])
            v[:, j] = shape * (v[:, j].sum() / max(shape.sum(), 1e-80))
        m = off + v @ spectra
        best = None
        for _ in range(3):
            for j in np.argsort(-v.sum(0)):
                use = spectra[j] > 0
                offset = np.maximum(
                    off[:, use], m[:, use] - v[:, j, None] * spectra[j, use]
                )
                v[:, j], _ = update_h(
                    t, y[:, use], offset, spectra[j, use], v[:, j], s.phi, price
                )
                m[:, use] = offset + v[:, j, None] * spectra[j, use]
            full = np.zeros_like(h0)
            full[:, keep] = v
            dd = cell_deviance_difference(y, m, old, s.phi).sum()
            dj = dd + price * (
                sum(face(t, col) for col in full.T) - df0 - last[full.sum(0) == 0].sum()
            )
            if best is None or dj < best[0]:
                best = dj, full.copy(), m.copy(), dd
        return best

    active = all_rows
    best = (0.0, h0, old, 0.0)
    for depth in range(len(candidates) + 1):
        subsets = [active] if depth == 0 else []
        subsets += [
            active[active != j] for j in candidates if j in active and len(active) > 1
        ]
        if not subsets:
            break
        proposals = ((u, refit(u, best[1])) for u in subsets)
        choice = min(
            (x for x in proposals if x[1] is not None),
            key=lambda x: x[1][0],
            default=None,
        )
        if choice is None:
            break
        keep, z = choice
        if z[0] >= best[0] - 1e-5:
            break
        best = z
        if len(keep) == len(active):
            break
        active = keep
    dj, v, m, dd = best
    if dj >= -1e-5:
        return None
    h[np.ix_(rows, ids)] = v
    mu[at] = m
    ns[ids[v.sum(0) == 0]] -= 1
    return dict(
        stage="sample overlap",
        source=" ".join(str(j + 1) for j in ids),
        receivers=" ".join(str(j + 1) for j in ids[v.sum(0) > 0]),
        GIC=dj,
        D=dd,
    )


def refine_overlaps(data, H, P, price):
    """Propose mixtures across samples, then compare joint profiles on each native scan grid."""
    if len(data) < 2 or len(P) < 2:
        return []
    mu = means(data, H, P)
    moves = _shared_exchanges(data, H, P, mu, price)
    count = np.array([h.sum(0) > 0 for h in H]).sum(0)
    for si, (s, h, m) in enumerate(zip(data, H, mu)):
        for group in _overlaps(h):
            result = _refit_window(s, h, P, m, group, count, price)
            if result is not None:
                moves.append(dict(result, sample=si + 1))
    return moves
