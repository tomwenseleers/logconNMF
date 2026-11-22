"""Three independently measured spectra contribute to the same mixed peak."""
from pathlib import Path
from types import SimpleNamespace
import sys
import numpy as np
import numba

sys.path.insert(0, sys.argv[1])
numba.config.CACHE_DIR = sys.argv[2]
from logconcounts.kernels import curve_error
from logconcounts.objective import means, score
from logconcounts.overlap_refit import refine_overlaps

t = np.arange(45, dtype=float)
shape = np.exp(-0.5 * ((t - 22) / 4) ** 2)
shape /= shape.sum()
pure = np.array([[0.7, 0.3, 0, 0, 0, 0], [0, 0, 0.6, 0.4, 0, 0], [0, 0, 0, 0, 0.8, 0.2]])
fraction = np.array([0.25, 0.35, 0.40])
p = np.vstack([pure, fraction @ pure])
samples, profiles = [], []
true_area = 250000 * fraction
for areas in (true_area, [80000, 0, 60000], [70000, 90000, 0], [0, 85000, 100000]):
    h = shape[:, None] * areas
    bg = np.broadcast_to(7.0, (len(t), p.shape[1])).copy()
    samples.append(SimpleNamespace(t=t, warp=t, y=bg + h @ pure, bg=bg, phi=2.0))
    profiles.append(np.column_stack([h, np.zeros_like(t)]))
profiles[0][:, :3] *= 0.12
profiles[0][:, 3] = 0.88 * true_area.sum() * shape
price = 15.0

def criterion():
    return score(samples, profiles, p, means(samples, profiles, p), price)["objective"]

before = criterion()
moves = refine_overlaps(samples, profiles, p, price)
assert criterion() < before
np.testing.assert_allclose(sum(m["GIC"] for m in moves), criterion() - before, atol=1e-6)
np.testing.assert_allclose(profiles[0].sum(0), np.r_[true_area, 0], rtol=0.005)
for s, h in zip(samples, profiles):
    assert all(curve_error(s.t, v) < 1e-8 for v in h.T)
np.testing.assert_allclose(p[:3].sum(1), np.ones(3))
assert refine_overlaps(samples, [h[:, :3].copy() for h in profiles], p[:3].copy(), price) == []
