"""Shared-spectrum factorization on native scan grids and original ion counts."""

from types import SimpleNamespace

import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits

from .anchors import anchor_dictionary
from .engine import (
    atoms,
    group_dictionary,
    lift,
    optimize,
    refine_support,
    spectral_block,
    trim,
)
from .kernels import curve_error, face, log_shape
from .local_spectra import local_spectra
from .objective import means, score
from .residual_pair import residual_pair


def _sample(part, groups, price):
    y = np.asarray(part["y"], dtype=float, order="C")
    z = np.array(part["background"], dtype=float, order="C")
    alpha = np.array(part["alpha"], dtype=float, order="C", copy=True)
    h = sparse.csc_matrix(part["H"])
    p = np.array(part["S"], dtype=float, order="C", copy=True)
    t = np.ravel(part["time"]).astype(float)
    warp = np.ravel(part["aligned"]).astype(float)
    grid = np.asarray(part["windows"], dtype=np.int64).reshape(-1, 3) - 1
    phi = float(part["dispersion"])
    if (
        y.shape != (len(t), len(groups))
        or h.shape != (len(t), len(p))
        or z.shape != (len(t), len(alpha))
        or alpha.shape[1] != p.shape[1]
        or len(grid) != len(p)
        or not np.isfinite(phi)
        or phi <= 0
        or np.any(np.diff(warp) <= 0)
    ):
        raise ValueError("Inconsistent count data, time coordinates or initial factors")
    s = SimpleNamespace(
        Y=y,
        Z=z,
        alpha=alpha,
        t=t,
        warp=warp,
        h=h,
        S=p,
        grid=grid,
        ap=grid[:, 2],
        phi=phi,
    )
    if part.get("grouped", False) or p.shape[1] != y.shape[1]:
        if p.shape[1] < max(groups) + 1:
            raise ValueError("Grouped spectrum does not match the channel partition")
        total = p.sum(1)
        count_h = h.toarray() * total
        original = lift([s], [count_h], p / np.maximum(total[:, None], 1e-100), groups)
        mu = means([s], [count_h], original)
        for _ in range(2):
            spectral_block([s], [count_h], original, mu, price)
        total = count_h.sum(0)
        keep = total > 0
        s.S = original[keep] * total[keep, None]
        s.h = sparse.csc_matrix(count_h[:, keep] / total[keep])
        s.grid = grid[keep]
        s.ap = s.grid[:, 2]
    else:
        s.y = y
        s.bg = np.maximum(z @ alpha + 1e-8, 1e-8)
    return s


def fit(parts, groups, price, cycles=12, exchanges=2, trace=False):
    """Return total-ion profiles, unit-sum spectra and the scored factor dimensions."""
    groups = np.asarray(groups, dtype=np.int64).reshape(-1) - 1
    price = float(price)
    with threadpool_limits(1):
        data = [_sample(x, groups, price) for x in parts]
        candidates = atoms(data)
        if not len(candidates.p):
            P = np.zeros((0, len(groups)))
            H = [np.zeros((len(s.t), 0)) for s in data]
            bands = np.zeros((0, 2))
        elif len(data) == 1:
            s = data[0]
            areas = s.S.sum(1)
            P = s.S / areas[:, None]
            H = [s.h.toarray() * areas]
            for k, h in enumerate(H[0].T):
                if curve_error(s.t, h) > 1e-9:
                    shape = log_shape(s.t, h)
                    H[0][:, k] = shape * (h.sum() / max(shape.sum(), 1e-100))
            bands = s.warp[s.grid[:, :2]]
            # A spectral match proposes a merge. Its actual count-model cost
            # decides whether distinct, closely spaced peaks stay separate.
            if len(groups) > 1:
                p, h, window = group_dictionary(data, candidates)
                if len(p) < len(P) and (
                    score(data, h, p, means(data, h, p), price)["objective"]
                    < score(data, H, P, means(data, H, P), price)["objective"]
                ):
                    P, H, bands = p, h, window
        else:
            P, H, bands = anchor_dictionary(data, candidates, price)
        del candidates
        history = [] if trace else None
        state = optimize(
            data,
            H,
            P,
            bands,
            price,
            int(cycles),
            history,
            "factor updates",
            active=True,
        )
        H, P, bands = trim(H, P, bands)
        if len(P):
            state = refine_support(
                data, H, P, bands, price, int(exchanges), history, state
            )
            before = state[1]["objective"]
            moves = local_spectra(data, H, P, bands, price)
            if len(data) > 1:
                P, bands, pairs = residual_pair(
                    data,
                    H,
                    P,
                    bands,
                    price,
                    replicated=True,
                    coexisting=True,
                    polished=True,
                )
                moves += pairs
        else:
            before, moves = state[1]["objective"], []
        H, P, bands = trim(H, P, bands)
        final = score(data, H, P, means(data, H, P), price)
        predicted = before + sum(z["GIC"] for z in moves)
        if not np.isclose(predicted, final["objective"], rtol=1e-8, atol=1e-4):
            raise ArithmeticError(
                "Conditional updates disagree with the original-ion objective"
            )
        if history is not None:
            history.append(
                dict(
                    stage="local joint refit",
                    cycle=1,
                    block="overlap",
                    channels=len(groups),
                    **final,
                )
            )
        # Return each identity in a stable elution order.
        order = np.argsort(bands.mean(axis=1), kind="stable")
        shape_df = []
        for s, h in zip(data, H):
            if any(curve_error(s.t, h[:, k]) > 1e-6 for k in order):
                raise ArithmeticError(
                    "Nonconcave log profile on native scan coordinates"
                )
            shape_df.append(np.asarray([face(s.t, h[:, k]) for k in order], dtype=int))
        return dict(
            P=P[order],
            H=[sparse.csc_matrix(h[:, order]) for h in H],
            alpha=[s.alpha for s in data],
            shape_df=shape_df,
            objective=final,
            trace=[] if history is None else history,
        )
