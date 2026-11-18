"""Recompute lattice energy errors from the saved fixed-grip solutions."""
import csv
from pathlib import Path
from statistics import mean


def verify(directory):
    with (directory / "results.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    fields = ("width", "notch", "patch_offset", "radius", "contrast")
    keys = {tuple(row[f] for f in fields) for row in rows}
    cases = [row for row in rows if float(row["contrast"]) != 1]
    if len(rows) != 336 or len(keys) != 336 or len(cases) != 288:
        raise ValueError("Expected 336 unique configurations, including 288 nonuniform cases")
    expected = {"local_control": "57.66", "uniform_transfer": "4.18", "tangent_transfer": "0.43"}
    for field, target in expected.items():
        errors = [100 * abs(float(r[field]) / float(r["released"]) - 1) for r in cases]
        observed = f"{mean(errors):.2f}"
        if observed != target:
            raise ValueError(f"{field} mean error was {observed}%, expected {target}%")
        if field == "tangent_transfer" and f"{max(errors):.2f}" != "5.96":
            raise ValueError("Maximum corrected error differs from the reported 5.96%")
    largest = max(float(row["identity_error"]) for row in rows)
    # Sparse solvers on different platforms can change cancellation at roundoff.
    if largest > 1e-10:
        raise ValueError(f"Energy identity residual too large: {largest}")
    print(f"Verified 336 configurations, all three error means and the worst corrected case. "
          f"Largest recorded energy identity residual: {largest:.3g}.")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
