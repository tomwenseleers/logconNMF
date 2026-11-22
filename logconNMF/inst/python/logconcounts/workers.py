"""One explicitly budgeted thread pool per fit. Workers never start another pool."""

from concurrent.futures import ThreadPoolExecutor
from contextvars import ContextVar
from contextlib import contextmanager

_context = ContextVar("joint_count_workers", default=(None, 1))


@contextmanager
def worker_pool(ncores):
    ncores = int(ncores)
    if ncores < 1:
        raise ValueError("ncores must be positive")
    pool = ThreadPoolExecutor(max_workers=ncores) if ncores > 1 else None
    token = _context.set((pool, ncores))
    try:
        yield
    finally:
        _context.reset(token)
        if pool is not None:
            pool.shutdown(wait=True)


def capacity():
    return _context.get()[1]


def ordered_map(fn, items):
    pool, _ = _context.get()
    if pool is None:
        return list(map(fn, items))

    def run(item):
        # Concurrent library calls and nested block updates get independent budgets.
        token = _context.set((None, 1))
        try:
            return fn(item)
        finally:
            _context.reset(token)

    return list(pool.map(run, items))


def initialize_samples(fn, parts, available):
    if not available or capacity() == 1:
        # Unknown memory headroom admits only one allocation-heavy lift at a time.
        return list(map(fn, parts))
    # Keep the native observations and baseline for every sample, plus each
    # active worker's count arrays, profile design and dense Gram workspace.
    retained = 0
    largest = 1
    for part in parts:
        n, m = part["y"].shape
        k = part["H"].shape[1]
        b = part["background"].shape[1]
        retained += 16 * n * m
        largest = max(largest, 8 * (8 * n * m + 4 * n * k + (k + b) ** 2))
    limit = max(1, min(capacity(), int((available * 0.75 - retained) // largest)))
    return [
        s
        for i in range(0, len(parts), limit)
        for s in ordered_map(fn, parts[i : i + limit])
    ]
