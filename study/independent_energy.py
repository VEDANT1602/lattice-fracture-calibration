"""Check the dipole observable against explicit energies before and after the tip cut."""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import argparse
import csv
from math import cos, exp, sin
from pathlib import Path
import numpy as np
from benchmark import Network


def stiffness(x, y, pattern):
    if pattern == "uniform":
        return 1.0
    if pattern == "wave":
        return exp(0.8*sin(0.45*x) + 0.4*cos(0.64*y))
    return 4.0 if (x - 0.65*y) > 1.3 else 0.25


def enumerate_bonds(n, pattern, release_tip):
    """Generate horizontal and vertical bonds, with the loaded tip at x=n/4."""
    bonds = []
    for y in range(n):
        for x in range(n):
            i = y*n + x
            if x+1 < n:
                bonds.append((i,i+1,stiffness(x+.5,y,pattern)))
            if y+1 == n:
                continue
            across = y == n//2-1
            if across and (x < n//4 or (release_tip and x == n//4)):
                continue
            k = 1.0 if across and x == n//4 else stiffness(x,y+.5,pattern)
            bonds.append((i,i+n,k))
    return bonds


def stored_energy(n, bonds):
    total_nodes = n*n
    full_stiffness = np.zeros((total_nodes, total_nodes))
    for i,j,k in bonds:
        full_stiffness[i,i] += k
        full_stiffness[j,j] += k
        full_stiffness[i,j] -= k
        full_stiffness[j,i] -= k
    u = np.zeros(total_nodes)
    u[-n:] = 1.0
    free = np.arange(n, total_nodes-n)
    A = full_stiffness[np.ix_(free,free)]
    u[free] = np.linalg.solve(A, -(full_stiffness@u)[free])
    energy = sum(0.5*k*(u[j]-u[i])**2 for i,j,k in bonds)
    return energy,u


def check():
    rows = []
    for n in (8,12,16):
        fast = Network(n)
        for pattern in ("uniform","wave","step"):
            pre,u = stored_energy(n, enumerate_bonds(n, pattern, False))
            post,_ = stored_energy(n, enumerate_bonds(n, pattern, True))
            tip_a = (n//2-1)*n+n//4
            tip_force = u[tip_a+n]-u[tip_a]
            measured = (pre-post)/tip_force**2
            edge_values = np.array([stiffness(x,y,pattern) for x,y in zip(fast.x,fast.y)])
            prediction = fast.exact(edge_values)
            difference = abs(measured/prediction-1)
            if difference > 1e-10:
                raise ValueError((n,pattern,measured,prediction,difference))
            rows.append(dict(n=n,pattern=pattern,energy_before=pre,energy_after=post,
                             pre_cut_force=tip_force,observed_normalized_release=measured,
                             dipole_prediction=prediction,relative_difference=difference))
    return rows


if __name__ == "__main__":
    opts = argparse.ArgumentParser()
    opts.add_argument("--out", type=Path,
                      default=Path(__file__).resolve().parent / "independent_energy.csv")
    path = opts.parse_args().out
    rows = check()
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("Nine energy checks; largest relative difference",max(r["relative_difference"] for r in rows))
