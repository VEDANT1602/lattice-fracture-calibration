"""Compare reconstructed outputs with the source repository's results.csv."""
import argparse
from csv import DictReader
from pathlib import Path

import numpy as np


def index(path, fields):
    with path.open() as stream:
        rows = list(DictReader(stream))
    keyed = {tuple(float(row[field]) for field in fields): row for row in rows}
    if len(rows) != 336 or len(keyed) != len(rows):
        raise ValueError(f"Unexpected count or duplicate cell in {path}")
    return keyed


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("archived_csv", type=Path)
args = parser.parse_args()
ours = index(Path(__file__).with_name("results.csv"), ("n", "notch", "offset", "radius", "contrast"))
archived = index(args.archived_csv, ("width", "notch", "patch_offset", "radius", "contrast"))
if ours.keys() != archived.keys():
    raise ValueError("Parameter grids differ")
for local, original in [("release", "released"), ("W", "local_control"),
                        ("grip", "critical_opening"), ("uniform_prediction", "uniform_transfer"),
                        ("corrected_prediction", "tangent_transfer")]:
    a = [float(ours[i][local]) for i in ours]
    b = [float(archived[i][original]) for i in ours]
    np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-13)
    print(f"{local}: maximum relative discrepancy {max(abs(np.asarray(a)/b-1)):.3g}")
