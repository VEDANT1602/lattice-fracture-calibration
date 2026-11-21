"""Check published tables against all per-field records without repeating timings."""
import csv
from math import isclose, sqrt
from pathlib import Path
from statistics import mean

HERE = Path(__file__).resolve().parent
FAMILIES = ("fixed_local", "distributed_mild", "distributed_strong")


def records(name):
    with (HERE / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def within(lower, value, upper):
    return lower - 1e-9 <= value <= upper + 1e-9


def main():
    setups = {int(row["n"]): row for row in records("setup.csv")}
    worst_hybrid = 0.0
    worst_woodbury = 0.0
    for n in (64, 208):
        timing = records(f"queries_{n}.csv")
        assert len(timing) == 144
        ref = float(setups[n]["uniform_output"])
        for family, expected_pass in zip(FAMILIES, (0, 48, 0), strict=True):
            rows = [r for r in timing if r["family"] == family]
            assert {int(r["seed"]) for r in rows} == set(range(81000, 81048))
            assert sum(int(r["screen_pass"]) for r in rows) == expected_pass
            for row in rows:
                exact = float(row["direct"])
                for method in ("screen", "tight_pcg", "cert_pcg", "tail_pcg"):
                    assert within(float(row[f"{method}_lower"]), exact, float(row[f"{method}_upper"]))
                for method in ("cert_pcg", "tail_pcg"):
                    lo, hi = (float(row[f"{method}_{side}"]) for side in ("lower", "upper"))
                    if method == "cert_pcg":
                        assert hi <= lo * 1.02**2 + 1e-11
                    else:
                        assert lo > ref or hi < ref
                error = abs(float(row["hybrid_value"]) / exact - 1)
                worst_hybrid = max(worst_hybrid, error)
                assert error <= 0.02
                if family == "fixed_local":
                    worst_woodbury = max(worst_woodbury, abs(float(row["woodbury"]) / exact - 1))
                else:
                    assert row["woodbury"] == ""
            if n == 208:
                print(f"{family}: average direct / tight PCG / 2% PCG / hybrid ms:",
                      *(f"{mean(float(r[col]) for r in rows):.3f}"
                        for col in ("direct_ms", "tight_pcg_ms", "cert_pcg_ms", "screen_hybrid_ms")))
    assert worst_woodbury <= 9e-16
    print(f"Max relative error, hybrid {100 * worst_hybrid:.5f}%, Woodbury {worst_woodbury:g}")

    draws = records("uq_draws.csv")
    summary = {r["family"]: r for r in records("uq_summary.csv")}
    assert len(draws) == 3 * 1024
    for family, count, fallback in zip(FAMILIES, (555, 508, 532), (0, 33, 170), strict=True):
        ref = float(summary[family]["threshold"])
        rows = [r for r in draws if r["family"] == family]
        assert {int(r["seed"]) for r in rows} == set(range(91000, 92024))
        assert sum(int(r["exceeds"]) for r in rows) == int(summary[family]["exceedances"]) == count
        assert sum(r["method"] == "pcg" for r in rows) == int(summary[family]["fallbacks"]) == fallback
        p = count / 1024
        assert isclose(float(summary[family]["mc_se"]), sqrt(p * (1 - p) / 1024))
        checked = 0
        for row in rows:
            lo, hi = float(row["lower"]), float(row["upper"])
            assert lo <= hi and (lo > ref or hi < ref)
            assert bool(int(row["exceeds"])) == (lo > ref)
            if row["exact_sentinel"]:
                exact = float(row["exact_sentinel"])
                checked += 1
                assert within(lo, exact, hi)
                assert (lo > ref) == (exact > ref)
        assert checked == int(summary[family]["independent_solver_checks"]) == 32
    print("Verified 288 full comparisons, 3,072 decisions and 96 independent solver checks.")


if __name__ == "__main__":
    main()
