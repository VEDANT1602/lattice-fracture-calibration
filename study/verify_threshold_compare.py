"""Recompute the reported threshold outcomes independently from saved per-field records."""
import csv
from collections import Counter
from math import isclose
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(name):
    with (ROOT/name).open(newline='') as f:
        return list(csv.DictReader(f))


def main():
    rows = read('threshold_fields.csv')
    sums = read('threshold_summary.csv')
    historical = {(r['family'],int(r['seed'])):r for r in read('uq_draws.csv')}
    keys = [(r['family'],int(r['seed'])) for r in rows]
    assert len(keys) == len(set(keys)) == 3072 and set(keys) == set(historical)
    sentinel_checks = 0
    for s in sums:
        family = s['family']
        cell = [r for r in rows if r['family']==family]
        assert len(cell)==int(s['fields'])==1024
        threshold = float(s['threshold'])
        expected = sum(float(r['exact_value'])>threshold for r in cell)
        assert expected==int(s['above'])
        paths = Counter(r['uniform_cholesky_path'] for r in cell)
        assert paths['pcg']==int(s['fallbacks'])
        for method in ('uniform_cholesky','uniform_superlu','field_cholesky'):
            times = [float(r[method+'_ms']) for r in cell]
            assert min(times)>0 and isclose(sum(times)/1000,float(s[method+'_seconds']), rel_tol=1e-12)
        for r in cell:
            value = float(r['exact_value'])
            old = historical[family,int(r['seed'])]
            if old['exact_sentinel']:
                assert isclose(float(old['exact_sentinel']),value,rel_tol=2e-11)
                sentinel_checks += 1
            for method in ('uniform_cholesky','uniform_superlu'):
                lo,hi = float(r[method+'_lower']),float(r[method+'_upper'])
                assert lo-1e-9<=value<=hi+1e-9 and (lo>threshold or hi<threshold)
                assert (lo>threshold)==(value>threshold)==bool(int(old['exceeds']))==bool(int(r['above']))
                assert r[method+'_path']==old['method']
    assert sentinel_checks==96 and len(sums)==3
    print('Verified both threshold pipelines against all 3,072 exact Cholesky solves, legacy classifications and 96 archived direct SuperLU checks.')


if __name__ == '__main__':
    main()
