"""Single-bond energy release with two independent equilibrium solves.

Requires Python 3, NumPy and SciPy. All quantities are dimensionless.
"""
import csv
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


ROOT = Path(__file__).resolve().parent


def geometry(n, a, center, radius):
    # n columns and n rows; a missing bonds on the horizontal midline.
    x, y = np.meshgrid(np.arange(n), np.arange(n))
    h = n // 2 - 1
    ix = np.arange(n * n).reshape(n, n)
    first = np.r_[ix[:, :-1].ravel(), ix[:-1].ravel()]
    second = np.r_[ix[:, 1:].ravel(), ix[1:].ravel()]
    keep = ~((second - first == n) & (y.flat[first] == h) & (x.flat[first] < a))
    first, second = first[keep], second[keep]
    tip = np.flatnonzero((first == ix[h, a]) & (second == ix[h + 1, a])).item()
    cx = (x.flat[first] + x.flat[second]) / 2
    cy = (y.flat[first] + y.flat[second]) / 2
    patch = (cx - a - center) ** 2 + (cy - h - 0.5) ** 2 <= radius ** 2
    rows = np.arange(len(first))
    incidence = sparse.coo_matrix((np.r_[-np.ones(len(first)), np.ones(len(first))],
                                  (np.r_[rows, rows], np.r_[first, second])),
                                  shape=(len(rows), n * n)).tocsr()
    boundary = (y.ravel() == 0) | (y.ravel() == n - 1)
    imposed = (y.ravel() == n - 1).astype(float)
    return incidence[:, ~boundary].tocsc(), incidence @ imposed, tip, patch


def release(geo, ratio):
    b, prescribed, tip, patch = geo
    k = np.where(patch, ratio, 1.)
    def equilibrium(weights):
        lu = splu((b.T @ b.multiply(weights[:, None])).tocsc())
        gains = lu.solve(-np.asarray(b.T @ (weights * prescribed)))
        ext = b @ gains + prescribed
        return lu, ext, 0.5 * np.dot(weights, ext ** 2)
    lu, ext, energy = equilibrium(k)
    h = b[tip:tip+1].toarray().ravel()
    c = h @ lu.solve(h)
    force = k[tip] * ext[tip]
    opening = 1 / abs(force)
    local = 0.5 / k[tip]
    # At the breaking load, compare the intact and re-equilibrated network.
    after = k.copy()
    after[tip] = 0
    _, _, e_after = equilibrium(after)
    released = (energy - e_after) * opening ** 2
    analytic = local / (1 - k[tip] * c)
    assert abs(released / analytic - 1) < 1e-8
    return released, local, analytic, opening


def main():
    fields = "width notch patch_offset radius contrast released local local_control uniform_transfer tangent_transfer identity_error critical_opening".split()
    with (ROOT / "results.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for n in (16, 24, 32, 48):
            for a in (n // 4, n // 2):
                for offset in (-4, 0, 4):
                    for radius in (1.25, 2.75):
                        geo = geometry(n, a, offset, radius)
                        def amp(r):
                            released, local, *_ = release(geo, r)
                            return np.log(released / local)
                        intercept = amp(1)
                        slope = (amp(np.exp(0.01)) - amp(np.exp(-0.01))) / 0.02
                        for r in (0.25, 0.5, 0.8, 1, 1.25, 2, 4):
                            truth, local, identity, opening = release(geo, r)
                            uniform = local * np.exp(intercept)
                            tangent = uniform * np.exp(slope * np.log(r))
                            w.writerow([n, a, offset, radius, r, truth, local, local, uniform, tangent, abs(identity / truth - 1), opening])
    records = np.genfromtxt(ROOT / "results.csv", delimiter=",", names=True)
    print("Contrast, local/uniform/tangent mean absolute relative errors (%):")
    predictions = ("local_control", "uniform_transfer", "tangent_transfer")
    for ratio in sorted(set(records["contrast"])):
        subset = records[records["contrast"] == ratio]
        errors = [100*np.mean(abs(subset[col]/subset["released"]-1)) for col in predictions]
        print(ratio, *np.round(errors, 3))
    print("Maximum relative residual of rank-one identity:", max(records["identity_error"]))


if __name__ == "__main__":
    main()
