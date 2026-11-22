"""Binary interchange with R. Matrices use R column order except raw observations."""

from pathlib import Path
import csv
import os
import sys


def available_memory():
    """Physical and commit headroom bound new Windows allocations independently."""
    if os.name == "nt":
        import ctypes
        class Status(ctypes.Structure):
            _fields_ = [("length", ctypes.c_uint32), ("load", ctypes.c_uint32)] + [
                (key, ctypes.c_uint64) for key in
                ("physical", "free", "commit", "uncommitted", "virtual", "address", "extended")
            ]
        status = Status()
        status.length = ctypes.sizeof(status)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            return min(status.free, status.uncommitted)
    else:
        try:
            return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        except (ValueError, OSError):
            pass
    return 0


def main(directory, cache):
    import numba

    numba.config.CACHE_DIR = str(cache)
    import numpy as np
    from scipy.sparse import csc_matrix
    from logconcounts.fit import fit

    def matrix(path, rows, cols):
        return np.fromfile(path, dtype="<f8").reshape(rows, cols, order="F")

    with (directory / "samples.csv").open(newline="") as f:
        meta = list(csv.DictReader(f))
    parts = []
    for i, row in enumerate(meta):
        d = directory / str(i + 1)
        n, m, b, k, q = [int(row[x]) for x in ("n", "m", "b", "k", "q")]
        y = np.memmap(row["observed"], dtype="<f8", mode="r", shape=(n, m), order="F")
        h = csc_matrix(
            (
                np.fromfile(d / "hx", dtype="<f8"),
                np.fromfile(d / "hi", dtype="<i4"),
                np.fromfile(d / "hp", dtype="<i4"),
            ),
            shape=(n, k),
        )
        parts.append(
            dict(
                y=y,
                grouped=row["grouped"] == "TRUE",
                H=h,
                S=matrix(d / "S", k, q),
                background=matrix(d / "B", n, b),
                alpha=matrix(d / "A", b, q),
                time=np.fromfile(d / "t", dtype="<f8"),
                aligned=np.fromfile(d / "warp", dtype="<f8"),
                dispersion=float(row["phi"]),
                windows=np.fromfile(d / "windows", dtype="<i4").reshape(
                    k, 3, order="F"
                ),
            )
        )
    groups = np.fromfile(directory / "groups", dtype="<i4")
    settings = meta[0]
    ans = fit(
        parts,
        groups,
        float(settings["price"]),
        int(settings["cycles"]),
        int(settings["exchanges"]),
        settings["trace"] == "TRUE",
    )
    ans["P"].T.astype("<f8").tofile(directory / "spectra")
    for i, (h, a, df) in enumerate(zip(ans["H"], ans["alpha"], ans["shape_df"])):
        d = directory / str(i + 1)
        h.indices.astype("<i4").tofile(d / "fi")
        h.indptr.astype("<i4").tofile(d / "fp")
        h.data.astype("<f8").tofile(d / "fx")
        a.T.astype("<f8").tofile(d / "fa")
        df.astype("<i4").tofile(d / "df")
    for name, rows in (("objective", [ans["objective"]]), ("trace", ans["trace"])):
        if rows:
            with (directory / (name + ".csv")).open("w", newline="") as f:
                out = csv.DictWriter(f, fieldnames=rows[0].keys())
                out.writeheader()
                out.writerows(rows)


if __name__ == "__main__":
    sys.dont_write_bytecode = True
    if sys.argv[1:] == ["--available-memory"]:
        print(available_memory())
        sys.exit(0)
    if sys.argv[1:] == ["--check"]:
        from importlib.util import find_spec

        sys.exit(
            int(
                not all(
                    find_spec(x) for x in ("numpy", "scipy", "numba", "threadpoolctl")
                )
            )
        )
    main(Path(sys.argv[1]), Path(sys.argv[2]))
