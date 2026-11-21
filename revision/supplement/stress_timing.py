"""Measure paired direct and hybrid runtimes for batches all known to require fallback."""
from time import perf_counter
from transfer import CutGraph, field, timed_batch, write_rows, np


def main():
    rows=[]
    for n in (16,32,48):
        for a in (n//4,n//2):
            g=CutGraph(n,a);g.calibrate()
            cases=[4.**(s*field(g,fam,seed)) for fam in ('iid_binary','smooth_binary')
                   for seed in range(12) for s in (-1,1)]
            assert all(g.fast_bounds(k)[1] > .02 for k in cases)
            rng=np.random.default_rng(77111+n+a)
            g.direct(cases[0])
            for rep in range(5):
                seq=[cases[i] for i in rng.permutation(len(cases))]
                methods={'direct':g.direct,'hybrid':lambda k:g.hybrid(k,.02)}
                row=dict(n=n, notch=a,rep=rep,cases=len(seq))
                for m in rng.permutation(list(methods)): row[m]=timed_batch(methods[m],seq)
                rows.append(row)
    write_rows('timings_high_contrast_ms.csv',rows)


if __name__=='__main__':
    main()
