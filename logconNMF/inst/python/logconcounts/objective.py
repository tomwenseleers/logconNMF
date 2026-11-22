import numpy as np
from scipy.special import xlogy
from .kernels import face


def means(data, H, P):
    return [np.maximum(s.bg + h @ P, 1e-8) for s, h in zip(data, H)]


def objective_terms(data, H, P, D, price):
    live = np.sum([h.sum(0) for h in H], axis=0) > 0
    spectral = (P[live] > 0).sum() - live.sum()
    faces = 0
    calls = 0
    for s, h in zip(data, H):
        ks = np.flatnonzero(h.sum(0) > 0)
        calls += len(ks)
        for k in ks:
            faces += face(s.t, h[:, k])
    return dict(
        D=float(D),
        objective=float(D + price * (faces + spectral)),
        rank=int(live.sum()),
        calls=calls,
        spectral_df=int(spectral),
        shape_df=int(faces),
    )


def score(data, H, P, mu, price):
    D = sum(2 / s.phi * np.sum(xlogy(s.y, s.y / m) - s.y + m) for s, m in zip(data, mu))
    return objective_terms(data, H, P, D, price)


def cell_deviance_difference(y, b, a, phi):
    # D(y,b) - D(y,a) without subtracting two large full objectives.
    d = b - a
    rel = d / a
    lr = np.log(b) - np.log(a)
    sel = np.abs(rel) < 0.5
    lr[sel] = np.log1p(rel[sel])
    return (2 / phi) * (d - y * lr)
