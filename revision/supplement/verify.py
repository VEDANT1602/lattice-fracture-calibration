"""Check model identities, summary values, and release order against the manuscript."""
from collections import Counter
from csv import DictReader
from pathlib import Path

import numpy as np

with Path(__file__).with_name("results.csv").open() as stream:
    rows = [{k: float(v) for k, v in row.items()} for row in DictReader(stream)]

# Rounded reported means for W, calibrated, and slope-corrected predictions.
contrast_means = {
    0.25: (55.73, 6.97, 1.69), 0.5: (56.62, 3.88, 0.37),
    0.8: (57.34, 1.31, 0.03), 1.0: (57.69, 0.0, 0.0),
    1.25: (58.04, 1.32, 0.02), 2.0: (58.73, 4.05, 0.15),
    4.0: (59.51, 7.57, 0.29),
}
errors = ("local_error_pct", "uniform_error_pct", "corrected_error_pct")
cells = {tuple(r[k] for k in ("n", "notch", "offset", "radius", "contrast")) for r in rows}
assert len(cells) == len(rows) == 336
for contrast, expected in contrast_means.items():
    group = [r for r in rows if r["contrast"] == contrast]
    assert len(group) == 48
    values = [np.mean([r[k] for r in group]) for k in errors]
    np.testing.assert_allclose(values, expected, rtol=0, atol=0.005)
assert max(r["identity_relative_residual"] for r in rows) < 1e-11
transfer = [r for r in rows if r["contrast"] != 1]
assert len(transfer) == 288
failing = [r for r in transfer if r["max_other_bond_force_at_tip_threshold"] > 1 + 1e-8]
assert Counter(r["contrast"] for r in failing) == {0.25: 16, 4: 7}
for k in errors:
    values = [r[k] for r in transfer]
    print(k, "(%, mean / median / P90 / max):", np.mean(values),
          *np.quantile(values, [.5, .9, 1], method="linear"))
violators = {id(r) for r in failing}
for label, subset in [("other bond first", failing), ("remaining", [r for r in transfer if id(r) not in violators])]:
    print(label, len(subset), [(k, np.mean([r[k] for r in subset])) for k in errors])
print("PASS: 336 configurations, 21 table entries, compliance identity and first-threshold check")
