"""Paired threshold timings using uniform LU or Cholesky and a field-specific full solve."""
import argparse
import copy
import csv
import json
import platform
from pathlib import Path

# Import before NumPy so the single-thread settings take effect.
from scaling import (FAMILIES, ROOT, SymbolicReuse, cholmod, cv_solve, cvxopt, np,
                     scipy, setup_screen, splu, timed)
from benchmark import Network, write_csv

METHODS = ('uniform_cholesky', 'uniform_superlu', 'field_cholesky')


def make_screen(base, uniform, backend):
    g = copy.copy(base)
    if backend == 'uniform_cholesky':
        setup_screen(g, uniform)
    else:
        g.base = splu(uniform)
        g.x0 = g.base.solve(g.h)
        g.r0 = float(g.h @ g.x0)
        g.setup_screen()
    g.setup_woodbury()
    return g


def classify(g, k, family):
    ref = (1 + g.r0) / 2
    if family == 'fixed_local':
        value = g.woodbury(k)
        return value, value, 'woodbury', 0
    lower, upper = g.screen(k, False)
    if lower > ref + 1e-11 or upper < ref - 1e-11:
        return lower, upper, 'screen', 0
    _, lower, upper, steps = g.pcg(k, 'tail')
    if lower <= ref <= upper:
        raise ArithmeticError('Decision not certified')
    return lower, upper, 'pcg', steps


def exact(g, k, sym):
    A = g.matrix(k)
    lower = sym.fill(A)
    cholmod.numeric(lower, sym.factor)
    x = cv_solve(sym.factor, g.h)
    return (1 + g.h @ x) / 2


def read_previous():
    with (ROOT/'uq_draws.csv').open() as stream:
        return {(r['family'], int(r['seed'])): r for r in csv.DictReader(stream)}


def run(out):
    out.mkdir(exist_ok=True, parents=True)
    n = 208
    cholmod.options['supernodal'] = 2
    # Exercise all three code paths on a small disposable network first.
    warm = Network(32)
    u = warm.matrix(warm.ones)
    for method in METHODS[:2]:
        small = make_screen(warm, u, method)
        for family in FAMILIES:
            classify(small, small.field(family, 91000), family)
    exact(warm, warm.ones, SymbolicReuse(u))
    g, geom_ms = timed(Network, n)
    u, assembly_ms = timed(g.matrix, g.ones)
    screens, setup = {}, {}
    for method in METHODS[:2]:
        screens[method], setup[method] = timed(make_screen, g, u, method)
    direct, setup['field_cholesky'] = timed(SymbolicReuse, u)
    ref = (1 + screens['uniform_cholesky'].r0) / 2
    previous = read_previous()
    rows = []
    rng = np.random.default_rng(712308)
    jobs = [(family, seed) for family in FAMILIES for seed in range(91000, 92024)]
    for i in rng.permutation(len(jobs)):
        family, seed = jobs[i]
        k = g.field(family, seed)
        row = dict(family=family, seed=seed)
        for method in rng.permutation(METHODS):
            if method == 'field_cholesky':
                value, elapsed = timed(exact, g, k, direct)
                row['exact_value'] = value
            else:
                (lo, hi, path, steps), elapsed = timed(classify, screens[method], k, family)
                row.update({method+'_lower': lo, method+'_upper': hi,
                            method+'_path': path, method+'_steps': steps})
            row[method+'_ms'] = elapsed
        value = row['exact_value']
        row['above'] = int(value > ref)
        old = previous[family, seed]
        assert int(old['exceeds']) == row['above']
        if old['exact_sentinel']:
            np.testing.assert_allclose(float(old['exact_sentinel']), value, rtol=2e-11)
        for m in METHODS[:2]:
            lo, hi = row[m+'_lower'], row[m+'_upper']
            assert lo-1e-9 <= value <= hi+1e-9 and (lo>ref or hi<ref)
            assert (lo>ref) == bool(row['above']) and row[m+'_path'] == old['method']
        rows.append(row)
        if len(rows) % 512 == 0:
            print(f'Verified {len(rows)} of {len(jobs)} fields', flush=True)
    cols = ['family','seed','above','exact_value','field_cholesky_ms']
    for method in METHODS[:2]:
        cols.extend(method+'_'+key for key in ('path','steps','lower','upper','ms'))
    rows.sort(key=lambda r:(r['family'],r['seed']))
    write_csv(out/'threshold_fields.csv', [{c:r[c] for c in cols} for r in rows])
    summary = []
    for family in FAMILIES:
        cell = [r for r in rows if r['family']==family]
        summary.append(dict(family=family,n=n,fields=len(cell),above=sum(r['above'] for r in cell),
                            fallbacks=sum(r['uniform_cholesky_path']=='pcg' for r in cell),
                            threshold=ref, **{m+'_seconds': sum(r[m+'_ms'] for r in cell)/1000 for m in METHODS}))
    write_csv(out/'threshold_summary.csv', summary)
    write_csv(out/'threshold_setup.csv',[dict(n=n,geometry_ms=geom_ms,uniform_assembly_ms=assembly_ms,
                                                 **{m+'_ms':v for m,v in setup.items()})])
    with open('/proc/cpuinfo') as f:
        cpu = next(line.split(':',1)[1].strip() for line in f if line.startswith('model name'))
    env = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
               cvxopt=cvxopt.__version__,platform=platform.platform(),cpu=cpu,
               randomization_seed=712308,cases=len(rows),single_thread=True)
    (out/'threshold_environment.json').write_text(json.dumps(env,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT)
    run(parser.parse_args().out)
