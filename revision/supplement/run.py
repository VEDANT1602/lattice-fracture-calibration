"""Reconstruct the fixed-grip lattice experiment from the supplied manuscript."""
from csv import DictWriter
from itertools import product
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


def lattice(n, notch, offset, radius):
    """Return bond incidence, patch membership and the tip bond index."""
    ids = np.arange(n * n).reshape(n, n)
    left = np.r_[ids[:, :-1].ravel(), ids[:-1, :].ravel()]
    right = np.r_[ids[:, 1:].ravel(), ids[1:, :].ravel()]
    midx = (left % n + right % n) / 2
    midy = (left // n + right // n) / 2
    midline = (n - 1) / 2
    keep = ~((midy == midline) & (midx < notch))
    left, right, midx, midy = [x[keep] for x in (left, right, midx, midy)]
    tip = int(np.flatnonzero((midy == midline) & (midx == notch)).item())
    patch = np.square(midx - notch - offset) + np.square(midy - midline) <= radius**2
    rows = np.arange(len(left))
    incidence = coo_matrix((np.r_[-np.ones_like(left), np.ones_like(left)],
                            (np.r_[rows, rows], np.r_[left, right])),
                           shape=(len(left), n * n)).tocsr()
    return incidence, patch, tip


def solve(n, geometry, contrast):
    bonds, patch, tip = geometry
    k = np.where(patch, contrast, 1.0)
    free = bonds[:, n:-n]
    imposed = np.asarray(bonds[:, -n:].sum(axis=1)).ravel()
    h = free.getrow(tip).toarray().ravel()

    def equilibrate(stiffness):
        weighted = free.T.multiply(stiffness).tocsr()
        lu = splu((weighted @ free).tocsc())
        return free @ lu.solve(-(weighted @ imposed)) + imposed, lu

    elongation, lu = equilibrate(k)
    tip_k = k[tip]
    grip = 1 / abs(tip_k * elongation[tip])
    w = 0.5 / tip_k
    c = h @ lu.solve(h)
    after = k.copy()
    after[tip] = 0
    cut_elongation, _ = equilibrate(after)
    release = 0.5 * grip**2 * (np.dot(k, elongation**2) - np.dot(after, cut_elongation**2))
    factor = 1 / (1 - tip_k * c)
    forces_at_cut = np.abs(k * elongation) * grip
    forces_at_cut[tip] = 0.0
    return grip, w, release, factor, float(np.max(forces_at_cut))


def experiment():
    for n in (16, 24, 32, 48):
        for notch, offset, radius in product((n // 4, n // 2), (-4, 0, 4), (1.25, 2.75)):
            geometry = lattice(n, notch, offset, radius)
            uniform_factor = solve(n, geometry, 1)[3]
            small = solve(n, geometry, np.exp(-0.01))[3]
            large = solve(n, geometry, np.exp(0.01))[3]
            slope = (np.log(large) - np.log(small)) / 0.02
            for contrast in (0.25, 0.5, 0.8, 1, 1.25, 2, 4):
                grip, w, release, factor, max_other_force = solve(n, geometry, contrast)
                pred_uniform = w * uniform_factor
                pred_corrected = pred_uniform * contrast ** slope
                yield dict(n=n, notch=notch, offset=offset, radius=radius, contrast=contrast,
                           grip=grip, W=w, release=release, rank_one_release=w*factor,
                           max_other_bond_force_at_tip_threshold=max_other_force,
                           identity_relative_residual=abs(w*factor/release-1),
                           uniform_prediction=pred_uniform, slope=slope,
                           corrected_prediction=pred_corrected,
                           local_error_pct=100*abs(w/release-1),
                           uniform_error_pct=100*abs(pred_uniform/release-1),
                           corrected_error_pct=100*abs(pred_corrected/release-1))


def main():
    rows = list(experiment())
    with Path(__file__).with_name('results.csv').open('w', newline='') as stream:
        writer = DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print('Reconstructed', len(rows), 'cells')
    for r in (0.25, 0.5, 0.8, 1, 1.25, 2, 4):
        print(r, [round(float(np.mean([x[k] for x in rows if x['contrast'] == r])), 2)
                  for k in ('local_error_pct', 'uniform_error_pct', 'corrected_error_pct')])
    transfer = [x for x in rows if x['contrast'] != 1]
    for k in ('local_error_pct', 'uniform_error_pct', 'corrected_error_pct'):
        errors = [x[k] for x in transfer]
        print(k, 'mean, median, 90th percentile, maximum:',
              np.mean(errors), *np.quantile(errors, [.5, .9, 1]))
    print('Maximum identity residual:', max(x['identity_relative_residual'] for x in rows))


if __name__ == '__main__':
    main()
