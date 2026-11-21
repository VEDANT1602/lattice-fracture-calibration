"""Independently validate the paired comparison and all table summaries."""
import csv
import json
from collections import defaultdict
from math import isclose
from pathlib import Path
from statistics import median

root = Path(__file__).resolve().parent
def read(name):
    with (root / name).open() as f:
        return list(csv.DictReader(f))

env = json.loads((root / "scaling_environment.json").read_text())
methods = {"chol_fresh", "chol_reuse", "superlu", "amg_fresh", "amg_reuse", "output_hybrid"}
groups = defaultdict(dict)
by_family = defaultdict(list)
for row in read("scaling_fields.csv"):
    n, sample, method = int(row["n"]), int(row["seed"]), row["method"]
    assert int(row["unknowns"]) == n*(n-2)
    key = n, row["family"], sample
    assert method not in groups[key]
    groups[key][method] = row
    by_family[n,row["family"],method].append(row)
    t = sum(float(row[s]) for s in ["assembly_ms","prep_ms","numeric_or_setup_ms","solve_ms"])
    assert isclose(t,float(row["total_ms"]),abs_tol=1.e-8)
    assert method=="output_hybrid" or float(row["relative_residual"]) <= 1.05e-10

delta, approximate = [], []
for group in groups.values():
    assert set(group) == methods
    ref = float(group["superlu"]["response"])
    for m in methods-{"output_hybrid"}:
        delta.append(abs(float(group[m]["response"])/ref-1))
    approximate.append(abs(float(group["output_hybrid"]["response"])/ref-1))
assert max(delta)<2.e-10 and max(approximate)<=.020000001
assert len(groups) == len(env["sizes"])*len(env["families"])*env["samples"]
seen = set()
for row in read("scaling_summary.csv"):
    key = (int(row["n"]),row["family"],row["method"])
    rr = by_family[key]
    assert len(rr) == int(row["fields"]) == env["samples"]
    for fld in ["total_ms","prep_ms","assembly_ms","numeric_or_setup_ms","solve_ms","iterations"]:
        assert isclose(float(row["median_"+fld]),median(float(s[fld]) for s in rr),abs_tol=1.e-8)
    assert key not in seen
    seen.add(key)
assert set(by_family) == seen
print(f"{len(groups)} fields verified; max point-solve difference {max(delta):.3g}; max 2% output error {max(approximate):.6%}.")
