"""One homogeneous dipole calibration, held-out transfer, and a timed exact fallback."""
import os
for _name in ["OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS"]:
    os.environ[_name] = "1"
import csv
import json
import platform
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
import scipy
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

OUT = Path(__file__).resolve().parent
FAMILIES = ("ellipse", "strip", "two_disks", "iid_binary", "smooth_binary")
AMPLITUDES = (1.25, 2.0, 4.0)


class CutGraph:
    """Fixed topology with sparse assembly indices shared by every tested method."""
    def __init__(self, n, notch):
        self.n, self.notch = n, notch
        ids = np.arange(n*n).reshape(n, n)
        i = np.r_[ids[:, :-1].ravel(), ids[:-1, :].ravel()]
        j = np.r_[ids[:, 1:].ravel(), ids[1:, :].ravel()]
        x = (i % n + j % n) / 2
        y = (i // n + j // n) / 2
        midline = (n-1)/2
        live = ~((y == midline) & (x < notch))
        i,j,x,y = [v[live] for v in [i,j,x,y]]
        self.i, self.j, self.x, self.y = i,j,x,y
        self.tip = int(np.flatnonzero((y == midline) & (x == notch)).item())
        # The cut network contains all intact springs except the prescribed spring.
        free = n*(n-2)
        self.free = free
        self.h = np.zeros(free)
        self.h[j[self.tip]-n] = 1
        self.h[i[self.tip]-n] = -1
        live_cut = np.arange(len(x)) != self.tip
        c, d = i-n, j-n
        rows, cols, refs, signs = [], [], [], []
        for first, second, sign in [(c,c,1), (d,d,1), (c,d,-1), (d,c,-1)]:
            keep = live_cut & (first>=0) & (first<free) & (second>=0) & (second<free)
            rows.extend(first[keep]); cols.extend(second[keep]); refs.extend(np.flatnonzero(keep))
            signs.extend(np.full(keep.sum(), sign))
        keys, inverse = np.unique(np.array(cols)*free+rows,return_inverse=True)
        self.indptr = np.r_[0, np.cumsum(np.bincount(keys//free,minlength=free))].astype(np.int32)
        self.indices = (keys%free).astype(np.int32)
        self.groups, self.refs, self.signs = inverse, np.asarray(refs), np.asarray(signs)

    def solve(self,k):
        vals = np.bincount(self.groups, weights=k[self.refs]*self.signs)
        K = csc_matrix((vals,self.indices,self.indptr),shape=(self.free,self.free))
        return splu(K).solve(self.h)

    def direct(self,k):
        u = self.solve(k)
        return 0.5*(1/k[self.tip] + self.h @ u)

    def calibrate(self):
        u = self.solve(np.ones(len(self.i)))
        self.r0 = self.h @ u
        full_u = np.pad(u,(self.n,self.n))
        voltage = full_u[self.j]-full_u[self.i]
        voltage[self.tip] = 0
        self.weight = voltage*voltage/self.r0
        np.testing.assert_allclose(self.weight.sum(),1,rtol=2e-12)

    def estimates(self,k):
        # q are uniform-state energy fractions of the alternative network.
        # The weighted geometric mean is the exponential first-order surrogate.
        base = 1/k[self.tip]
        lower = 0.5*(base+self.r0/(self.weight @ k))
        upper = 0.5*(base+self.r0*(self.weight @ (1/k)))
        center = np.sqrt(lower*upper)
        first = 0.5*(base+self.r0*np.exp(-self.weight @ np.log(k)))
        return (0.5*(base+self.r0),first,center,lower,upper)

    def fast_bounds(self,k):
        base = 1/k[self.tip]
        lower = 0.5*(base+self.r0/(self.weight @ k))
        upper = 0.5*(base+self.r0*(self.weight @ (1/k)))
        return np.sqrt(lower*upper), np.sqrt(upper/lower)-1

    def hybrid(self,k,tolerance):
        estimate,bound = self.fast_bounds(k)
        return estimate if bound <= tolerance else self.direct(k)


def field(g,family,seed):
    rng=np.random.default_rng(307498+seed)
    dx,dy = g.x-g.notch,g.y-(g.n-1)/2
    # Shared draw across geometries for patch shape parameters; independent keyed grids for noise.
    cx,cy = rng.uniform(-3,3,2)
    theta = rng.uniform(0,np.pi)
    c,s=np.cos(theta),np.sin(theta)
    u = c*(dx-cx)+s*(dy-cy)
    v = -s*(dx-cx)+c*(dy-cy)
    if family == "ellipse":
        return (((u/3.1)**2+(v/1.45)**2)<=1).astype(float)
    if family == "strip":
        return ((abs(u)<=5) & (abs(v)<=0.85)).astype(float)
    if family == "two_disks":
        return (((u-2.4)**2+v**2<=1.4**2) | ((u+2.4)**2+v**2<=1.8**2)).astype(float)
    if family == "iid_binary":
        return rng.choice([-1.,1.],len(g.x))
    nodes = gaussian_filter(rng.standard_normal((g.n,g.n)),sigma=2.,mode="reflect")
    sample = map_coordinates(nodes,[g.y,g.x],order=1,mode="nearest")
    return np.where(sample>=np.median(sample),1.,-1.)


def records(g):
    for family in FAMILIES:
        for seed in range(12):
            z=field(g,family,seed)
            for ratio in AMPLITUDES:
                for direction in (-1,1):
                    yield (family,seed,ratio,direction), np.power(ratio,direction*z)


def timed_batch(fn,inputs):
    start=perf_counter()
    checksum = sum(fn(k) for k in inputs)
    elapsed=perf_counter()-start
    if not np.isfinite(checksum): raise ArithmeticError
    return elapsed*1000/len(inputs)


def write_rows(name,rows):
    with (OUT/name).open("w",newline="") as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    rows,timings, checks = [],[],[]
    for n in (16,32,48):
        for a in (n//4,n//2):
            g=CutGraph(n,a)
            t0=perf_counter();g.calibrate();first_setup=1000*(perf_counter()-t0)
            setup=[]
            for _ in range(9):
                t0=perf_counter();g.calibrate();setup.append(1000*(perf_counter()-t0))
            rng=np.random.default_rng(8800+n+a)
            # Keep the 360 inputs in memory for identical workloads in each timing block.
            cases=list(records(g))
            for (family,seed,ratio,sign),k in cases:
                reference=g.direct(k)
                base,linear,center,lo,hi=g.estimates(k)
                rows.append(dict(n=n,notch=a,family=family,seed=seed,ratio=ratio,sign=sign,
                           exact=reference,uniform=base,log_sensitivity=linear,center=center,lower=lo,upper=hi,
                           tip_k=k[g.tip],guaranteed_relative_bound=np.sqrt(hi/lo)-1))
                if seed == 0 and ratio == 4 and sign == 1:
                    # Independent before/after energy calculation by incidence assembly.
                    from scipy.sparse import coo_matrix
                    edge=np.arange(len(k))
                    B=coo_matrix((np.r_[-np.ones(len(k)),np.ones(len(k))],
                          (np.r_[edge,edge],np.r_[g.i,g.j])),shape=(len(k),n*n)).tocsr()
                    free=B[:,n:-n]
                    fixed=np.asarray(B[:,-n:].sum(axis=1)).ravel()
                    def energy(ks):
                        weighted=free.T.multiply(ks).tocsr()
                        u=splu((weighted@free).tocsc()).solve(-(weighted@fixed))
                        e=free@u+fixed
                        return .5*np.dot(ks,e*e),ks[g.tip]*e[g.tip]
                    pre,force=energy(k); ks=k.copy();ks[g.tip]=0
                    post,_=energy(ks)
                    difference=(pre-post)/force**2
                    checks.append(dict(n=n,notch=a,family=family,
                                       relative_difference=abs(difference/reference-1)))
            # Sample 90 of 360 cases without replacement; reshuffle these 90 in each repetition.
            # No array generation or randomization is inside the timed blocks.
            selected=rng.permutation(len(cases))[:90]
            for rep in range(5):
                order=rng.permutation(selected)
                inputs=[cases[j][1] for j in order]
                methods={"direct":g.direct, "bounds_only":lambda k:g.fast_bounds(k)[0],
                         "hybrid_2pct":lambda k:g.hybrid(k,.02)}
                row=dict(n=n,notch=a,rep=rep,cases=len(inputs),setup_first_ms=first_setup,
                         setup_median_ms=np.median(setup))
                for method in rng.permutation(list(methods)):
                    row[method]=timed_batch(methods[method],inputs)
                timings.append(row)
    write_rows("held_out.csv",rows)
    # Randomized timing order changes key insertion order, CSV fields use first row's field order.
    write_rows("timings_ms.csv",timings);write_rows("energy_checks.csv",checks)
    (OUT/"environment.json").write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,
        scipy=scipy.__version__,processor=platform.processor(),platform=platform.platform(),
        thread_limit=1),indent=2)+"\n")


if __name__ == "__main__":
    main()
