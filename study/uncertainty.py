"""Probability of exceeding the uniform release on fresh random-field draws."""
import csv
import json
from pathlib import Path
from time import perf_counter
from benchmark import Network,FAMILIES,ROOT,write_csv,np,scipy

N=208
SIZE=1024

def main():
    started=perf_counter()
    g=Network(N)
    g.factor_base();g.setup_screen();g.setup_woodbury()
    setup=perf_counter()-started
    reference=(1+g.r0)/2
    rng=np.random.default_rng(839017)
    sentinel=set(rng.choice(SIZE,32,replace=False))
    rows=[]; summaries=[]
    for family in FAMILIES:
        tgen=tsolve=tcheck=0
        successes=escalated=checked=0
        for i in range(SIZE):
            seed=91000+i
            t=perf_counter();k=g.field(family,seed);tgen+=perf_counter()-t
            t=perf_counter()
            if family=='fixed_local':
                result=g.woodbury(k)
                lo=hi=result
                path='woodbury'
            else:
                lo,hi=g.screen(k,False)
                if lo<=reference+1e-11 and hi>=reference-1e-11:
                    _,lo,hi,step=g.pcg(k,'tail')
                    path='pcg';escalated+=1
                else:path='screen'
            if lo>reference:exceeds=1
            elif hi<reference:exceeds=0
            else:raise RuntimeError('undecided threshold; run exact solver')
            tsolve+=perf_counter()-t
            exact=''
            if i in sentinel:
                t=perf_counter()
                exact=g.exact(k)
                assert (exact>reference)==bool(exceeds)
                assert lo-1e-9<=exact<=hi+1e-9
                checked+=1;tcheck+=perf_counter()-t
            successes+=exceeds
            rows.append(dict(family=family,seed=seed,lower=lo,upper=hi,method=path,exceeds=exceeds,exact_sentinel=exact))
        p=successes/SIZE
        summaries.append(dict(family=family,n=N,replicates=SIZE,threshold=reference,exceedances=successes,probability=p,
            mc_se=np.sqrt(p*(1-p)/SIZE),fallbacks=escalated,independent_solver_checks=checked,common_setup_s=setup,
            generation_s=tgen,query_s=tsolve,sentinel_direct_s=tcheck))
        print(family,successes,escalated,flush=True)
    write_csv(ROOT/'uq_draws.csv',rows)
    write_csv(ROOT/'uq_summary.csv',summaries)
if __name__=='__main__':main()
