"""Check the archived transfer predictions, reference solves and timing observations."""
import csv
import math
from collections import Counter
from pathlib import Path
from statistics import mean, median

ROOT = Path(__file__).resolve().parent

def read(name):
    with (ROOT/name).open(newline="") as stream:
        return list(csv.DictReader(stream))


rows=read("held_out.csv")
keys=("n","notch","family","seed","ratio","sign")
assert len(rows)==len({tuple(r[k] for k in keys) for r in rows})==2160
assert set(Counter(r["family"] for r in rows).values())=={432}
summary=[]
for tol in (.005,.01,.02,.05,.1):
    accepted=[]
    for r in rows:
        lo,ex,hi,proxy=(float(r[k]) for k in ("lower","exact","upper","center"))
        assert lo <= ex*(1+1e-10) and hi >= ex*(1-1e-10)
        assert math.isclose(proxy,math.sqrt(lo*hi),rel_tol=1e-12)
        cert=math.sqrt(hi/lo)-1
        assert math.isclose(cert,float(r["guaranteed_relative_bound"]),abs_tol=1e-12)
        if cert<=tol:
            error=abs(proxy/ex-1)
            assert error<=tol+1e-12
            accepted.append(error)
    summary.append((tol,len(accepted),100*max(accepted),100*sum(accepted)/len(rows)))
assert [r[1] for r in summary]==[564,772,1037,1443,1700]
print("tolerance / direct solves avoided / worst hybrid error (%) / mean hybrid error (%)")
for line in summary: print(*line)
for r in read('energy_checks.csv'):
    assert float(r['relative_difference'])<1e-10
assert len(read('energy_checks.csv'))==30
print("All 30 two-equilibrium reference checks agree.")
t=read("timings_ms.csv"); h=read("timings_high_contrast_ms.csv")
assert len(t)==len(h)==30
print("n, notch, direct / hybrid [ms], setup [ms], speedup, amortization [queries]")
for n,a in sorted({(int(r['n']),int(r['notch'])) for r in t}):
    group=[r for r in t if int(r['n'])==n and int(r['notch'])==a]
    assert len(group)==5
    get=lambda key: median([float(r[key]) for r in group])
    direct,hybrid,setup=[get(k) for k in ('direct','hybrid_2pct','setup_median_ms')]
    assert min(direct,hybrid,setup)>0
    print(n,a,round(direct,3),round(hybrid,3),round(setup,3),
          round(direct/hybrid,2),math.ceil(setup/(direct-hybrid)) if direct>hybrid else 'no break-even')
