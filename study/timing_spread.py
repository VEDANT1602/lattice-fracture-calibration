"""Quartiles and setup amortization from the published five-size timing panel."""
import csv
from collections import defaultdict
from math import ceil
from pathlib import Path
from statistics import mean, median, quantiles

ROOT = Path(__file__).resolve().parent


def read(name):
    with (ROOT/name).open(newline='') as f:
        return list(csv.DictReader(f))


def write(name, rows):
    with (ROOT/name).open('w',newline='') as f:
        out = csv.DictWriter(f,fieldnames=rows[0])
        out.writeheader()
        out.writerows(rows)


def main():
    cells = defaultdict(dict)
    for r in read('scaling_fields.csv'):
        key = int(r['n']),r['family'],r['method']
        seed = int(r['seed'])
        assert seed not in cells[key]
        cells[key][seed] = float(r['total_ms'])
    spread = []
    for (n,family,method), field_times in sorted(cells.items()):
        assert sorted(field_times) == list(range(183000,183008))
        times = list(field_times.values())
        q1, q2, q3 = quantiles(times, method='inclusive')
        spread.append(dict(n=n,family=family,method=method,fields=len(times),
                           q25_ms=q1,median_ms=q2,q75_ms=q3,iqr_width_ms=q3-q1))
    write('scaling_quartiles.csv',spread)
    amortization = []
    for row in read('scaling_setup.csv'):
        n=int(row['n'])
        premium = float(row['uniform_screen_and_factor_ms']) - float(row['reusable_symbolic_ms'])
        for family in ('fixed_local','distributed_mild','distributed_strong'):
            direct = cells[n,family,'chol_reuse']
            scalar = cells[n,family,'output_hybrid']
            gap = median(direct.values()) - median(scalar.values())
            gains = [direct[i]-scalar[i] for i in sorted(direct)]
            first = next((i for i in range(1,len(gains)+1) if sum(gains[:i]) >= premium), '')
            amortization.append(dict(n=n,family=family,extra_setup_ms=premium,
                                     difference_of_medians_ms=gap,
                                     queries_at_medians=ceil(premium/gap) if gap>0 else '',
                                     mean_paired_saving_ms=mean(gains),
                                     first_crossing_in_ascending_seed_order=first,
                                     total_saving_net_of_setup_for_eight_queries_ms=sum(gains)-premium))
    write('scaling_amortization.csv',amortization)


if __name__ == '__main__':
    main()
