"""Plot output-target and direct-solver measurements for strong distributed fields."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
with (root / "scaling_summary.csv").open() as stream:
    rows = [r for r in csv.DictReader(stream) if r["family"] == "distributed_strong"]
fig, ax = plt.subplots(figsize=(7.2, 3.7))
styles = [
    ("amg_reuse", "Uniform AMG + PCG", "#79528c", "P"),
    ("amg_fresh", "Field AMG + PCG", "#275a92", "v"),
    ("superlu", "Fresh SuperLU", "#656565", "^"),
    ("chol_fresh", "Fresh CHOLMOD", "#168078", "o"),
    ("chol_reuse", "Reused CHOLMOD pattern", "#ba611f", "s"),
    ("output_hybrid", "2% scalar certificate", "#232733", "D"),
]
for name, label, color, marker in styles:
    records = sorted((r for r in rows if r["method"] == name), key=lambda r: int(r["n"]))
    ax.plot([int(r["n"]) for r in records], [float(r["median_total_ms"]) for r in records],
            color=color, linestyle="--" if name=="output_hybrid" else "-", marker=marker,
            markersize=4, linewidth=1.6, label=label)
ax.set_xlabel("Lattice side length (nodes)")
ax.set_ylabel("Median query time (ms)")
ax.set_yscale("log")
ax.set_xticks(sorted({int(r["n"]) for r in rows}))
ax.spines[["top","right"]].set_visible(False)
ax.grid(axis="y", alpha=.23)
ax.legend(loc="upper left", fontsize=7, ncols=2, frameon=False)
fig.tight_layout()
fig.savefig(root / "scaling.png", dpi=220)
