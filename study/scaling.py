"""Compare direct and multigrid methods on identical finite-network fields."""
import os
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
import json
import platform
import sys
from pathlib import Path
from statistics import median
from time import perf_counter_ns

import cvxopt
from cvxopt import cholmod
import numpy as np
import pyamg
import scipy
from scipy.sparse import tril
from scipy.sparse.linalg import cg, splu
from benchmark import Network, write_csv

ROOT = Path(__file__).resolve().parent
FAMILIES = ("fixed_local", "distributed_mild", "distributed_strong")
METHODS = ("superlu", "chol_fresh", "chol_reuse", "amg_fresh", "amg_reuse")
RTOL = 1e-10


def timed(fn, *args):
    start = perf_counter_ns()
    x = fn(*args)
    return x, (perf_counter_ns() - start) / 1e6


def lower_triangle(matrix):
    low = tril(matrix, format="coo")
    return cvxopt.spmatrix(low.data, low.row, low.col, size=matrix.shape, tc="d")


class SymbolicReuse:
    """Reuse symbolic ordering and the sparse lower-triangle storage."""
    def __init__(self, A):
        col = np.repeat(np.arange(A.shape[1]), np.diff(A.indptr))
        self.keep = np.flatnonzero(A.indices >= col)
        self.lower = lower_triangle(A)
        np.testing.assert_array_equal(np.array(self.lower.V).ravel(), A.data[self.keep])
        self.factor = cholmod.symbolic(self.lower)

    def fill(self, A):
        self.lower.V = cvxopt.matrix(A.data[self.keep])
        return self.lower


class UniformFactor:
    def __init__(self, A):
        lo = lower_triangle(A)
        self.factor = cholmod.symbolic(lo)
        cholmod.numeric(lo, self.factor)

    def solve(self, b):
        return cv_solve(self.factor, b)


def setup_screen(g, uniform):
    g.base = UniformFactor(uniform)
    g.x0 = g.base.solve(g.h)
    g.r0 = g.h @ g.x0
    g.setup_screen()


def output_query(g, k, local):
    lo, hi = g.screen(k, local)
    if hi <= lo*1.02**2:
        return np.sqrt(lo*hi), 0
    value, _, _, steps = g.pcg(k, "output")
    return value, steps


def cv_solve(factor, h):
    b = cvxopt.matrix(h)
    cholmod.solve(factor, b)
    return np.array(b).ravel()


def hierarchy(matrix):
    return pyamg.ruge_stuben_solver(matrix, max_coarse=10, coarse_solver="pinv")


def iterative(matrix, ml, h):
    visits = []
    vcycle = ml.aspreconditioner(cycle="V")
    x, info = cg(matrix, h, M=vcycle, rtol=RTOL, atol=0, maxiter=350,
                 callback=lambda _: visits.append(None))
    if info:
        raise ArithmeticError(f"PCG stopped with info={info}")
    return x, len(visits)


def csr_bytes(A):
    return A.data.nbytes + A.indices.nbytes + A.indptr.nbytes


def amg_bytes(ml):
    """Sparse arrays in the hierarchy, counting each allocated backing array once."""
    buffers = {}
    for level in ml.levels:
        for name in ("A", "P", "R"):
            a = getattr(level, name, None)
            if a is None:
                continue
            for arr in (a.data, a.indices, a.indptr):
                backing = arr
                while isinstance(backing.base, np.ndarray):
                    backing = backing.base
                buffers[id(backing)] = backing.nbytes
    return sum(buffers.values())


def memory_record(A):
    def sparse_chol_factor():
        lower = lower_triangle(A)
        f = cholmod.symbolic(lower)
        cholmod.numeric(lower, f)
        return cholmod.getfactor(f)
    L = sparse_chol_factor()
    ptr, idx, val = L.CCS
    chol_bytes = sum(np.asarray(x).nbytes for x in (ptr, idx, val))
    lu = splu(A)
    multilevel = hierarchy(A.tocsr())
    return dict(matrix_csc_bytes=csr_bytes(A), chol_L_nnz=len(L),
                chol_exported_L_bytes=chol_bytes, lu_LU_nnz=lu.L.nnz+lu.U.nnz,
                lu_exported_LU_bytes=csr_bytes(lu.L)+csr_bytes(lu.U),
                amg_levels=len(multilevel.levels), amg_APR_bytes=amg_bytes(multilevel))


def solve_field(method, A, h, sym, uniform_amg):
    prep = numeric = 0.0
    iterations = 0
    if method == "superlu":
        lu, numeric = timed(splu, A)
        x, solve = timed(lu.solve, h)
    elif method.startswith("chol"):
        if method == "chol_fresh":
            low, prep = timed(lower_triangle, A)
            f, symbolic = timed(cholmod.symbolic, low)
            prep += symbolic
        else:
            low, prep = timed(sym.fill, A)
            f = sym.factor
        _, numeric = timed(cholmod.numeric, low, f)
        x, solve = timed(cv_solve, f, h)
    else:
        csr, prep = timed(A.tocsr)
        ml = uniform_amg
        if method == "amg_fresh":
            ml, numeric = timed(hierarchy, csr)
        (x, iterations), solve = timed(iterative, csr, ml, h)
    relres = np.linalg.norm(A@x-h) / np.linalg.norm(h)
    if relres > 1.05*RTOL:
        raise ArithmeticError(f"{method}: relative residual {relres}")
    return dict(prep_ms=prep, numeric_or_setup_ms=numeric, solve_ms=solve,
                iterations=iterations, relative_residual=float(relres),
                response=float((1+h@x)/2))


def sample(n, count):
    g, geom_time = timed(Network, n)
    uniform, uniform_assembly = timed(g.matrix, g.ones)
    symbolic, symbolic_time = timed(SymbolicReuse, uniform)
    _, screen_time = timed(setup_screen, g, uniform)
    ml, amg_time = timed(lambda: hierarchy(uniform.tocsr()))
    # Trigger the lazy coarse-level pseudoinverse once, charging it to the reusable setup.
    _, coarse_warmup = timed(lambda: ml.solve(g.h, maxiter=1, cycle="V", tol=0))
    setup = dict(n=n, unknowns=g.dim, geometry_ms=geom_time,
                 uniform_assembly_ms=uniform_assembly, reusable_symbolic_ms=symbolic_time,
                 uniform_amg_ms=amg_time+coarse_warmup, uniform_screen_and_factor_ms=screen_time)
    rng = np.random.default_rng(247135+n)
    records = []
    memory = []
    jobs = [(family, 183000+j) for family in FAMILIES for j in range(count)]
    for order in rng.permutation(len(jobs)):
        family, seed = jobs[order]
        k = g.field(family, seed)
        A, assembly = timed(g.matrix, k)
        values = {}
        for method in rng.permutation((*METHODS, "output_hybrid")):
            if method == "output_hybrid":
                (value, steps), elapsed = timed(output_query, g, k, family == "fixed_local")
                # The directly timed pipeline assembles internally and only on fallback.
                s = dict(prep_ms=0., numeric_or_setup_ms=0., solve_ms=elapsed,
                         iterations=steps, relative_residual="", response=float(value),
                         total_ms=elapsed, n=n, unknowns=g.dim, family=family, seed=seed,
                         method=str(method), assembly_ms=0.)
                records.append(s)
                continue
            s = solve_field(method, A, g.h, symbolic, ml)
            s["total_ms"] = assembly + s["prep_ms"] + s["numeric_or_setup_ms"] + s["solve_ms"]
            s.update(n=n, unknowns=g.dim, family=family, seed=seed,
                     method=str(method), assembly_ms=assembly)
            values[method] = s["response"]
            records.append(s)
        ref = values["superlu"]
        if abs(value/ref-1) > 0.020000001:
            raise ArithmeticError("Scalar certificate failed")
        if max(abs(v/ref-1) for v in values.values()) > 2e-10:
            raise ArithmeticError(f"Inconsistent response: {values}")
        if seed == 183000:
            memory.append(dict(n=n, family=family, **memory_record(A)))
    return setup, records, memory


def summarize(records):
    import itertools
    out = []
    for key, rows in itertools.groupby(sorted(records,key=lambda r:(r["n"],r["family"],r["method"])),
                                       key=lambda r:(r["n"],r["family"],r["method"])):
        a = list(rows)
        d = dict(n=key[0], family=key[1], method=key[2], fields=len(a))
        for field in ("total_ms", "assembly_ms", "prep_ms", "numeric_or_setup_ms", "solve_ms", "iterations"):
            d["median_"+field] = median(r[field] for r in a)
        d["min_total_ms"] = min(r["total_ms"] for r in a)
        d["max_total_ms"] = max(r["total_ms"] for r in a)
        out.append(d)
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--sizes", type=int, nargs="+", default=[64, 96, 144, 208, 288])
    cli.add_argument("--count", type=int, default=8)
    cli.add_argument("--out", type=Path, default=ROOT)
    opts = cli.parse_args()
    opts.out.mkdir(parents=True, exist_ok=True)
    # Explicit supernodal LL' (default ordering), lower triangle only.
    cholmod.options["supernodal"] = 2
    # Warm imported paths before the first measured specimen.
    sample(24, 1)
    setups, records, memory = [], [], []
    for n in opts.sizes:
        s, r, m = sample(n, opts.count)
        setups.append(s)
        records.extend(r)
        memory.extend(m)
        for name, rows in [("scaling_fields", records),("scaling_setup", setups),
                           ("scaling_memory", memory),("scaling_summary", summarize(records))]:
            write_csv(opts.out/(name+".csv"), rows)
        print(f"Completed n={n}", flush=True)
    env = dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
               cvxopt=cvxopt.__version__, pyamg=pyamg.__version__, platform=platform.platform(),
               single_thread=True, families=list(FAMILIES), sizes=opts.sizes, samples=opts.count)
    (opts.out/"scaling_environment.json").write_text(json.dumps(env, indent=2)+"\n")
