# Repeated stiffness queries on a cut spring network

This folder contains the scripts, input model and results for the scalar network benchmark. All fields are generated from seeds. No external material measurements are used.

## Run

Use CPython 3.12 on Linux x86_64. The locked binary wheels reproduce the NumPy and SciPy versions recorded in environment.json.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install --require-hashes -r requirements.lock
python benchmark.py
python uncertainty.py
python verify.py
```

The first two commands that invoke study scripts regenerate CSV files in this directory. Make a separate copy first if preserving the distributed timing records. To check the archived tables without rerunning timing measurements, use only `python verify.py`. That command uses the standard library. Execution times vary with the hardware, load, numerical libraries and their build. The studies request one thread each for OpenBLAS, MKL and OpenMP.

## Data and estimands

* `setup.csv` contains geometry setup shared by all approaches, uniform factorization and right-hand-side solve required by all recycled approaches, extra trial-field setup and extra 16-column Woodbury setup. Direct SuperLU requires geometry only.
* `queries_64.csv` and `queries_208.csv` contain 48 fresh fields per family and size. Solver timings include field-dependent matrix assembly for those that need it. Field generation is outside each compared timer. The fixed-support shortcuts exploit the 16 prescribed edge indices. The dense-rank Woodbury storage estimate is the size of one dense edge-by-edge matrix, not a lower bound on every possible implementation. Direct LU solves the sparse system from scratch each time. The 2% PCG implementation and tight PCG use the same preconditioner and initial guess. The screened hybrid's point-query column sums individual per-field component times, including fallback when required. Independent runs should reproduce numerical conclusions within double-precision tolerance. Latencies and bit patterns can differ.
* `uq_draws.csv` gives the certified interval, decision and optional direct-solver check for every field in the many-query experiment. `uq_summary.csv` separates per-family input generation, actual query pipeline and sentinel-check times from shared amortized setup. Use the threshold decision directly when estimating an exceedance probability. A certified 2% point estimate alone can give the wrong indicator.
* `CITATION-CHECKS.md` identifies the source sections underlying the paper. Timing tests compare sparse LU, known-support Woodbury, uniform-factor PCG and two forms of stopping. Multigrid and continuum approximations have no timing entries in these experiments.

## Stipulated design choices

An even square side length of 64 provides a smaller comparison; 208 provides 42,848 free nodes and 86,059 surviving bonds. Bond spacing, intact tip stiffness and baseline stiffness are one. Both grip rows have prescribed displacement. Sixteen nearest surviving bonds define the known-support case, keeping the modification low rank without assuming all perturbations are weak. Their stiffnesses vary independently and log-uniformly between a quarter and four times unit stiffness. Distributed fields threshold Gaussian noise smoothed over three lattice spacings, producing connected clusters. Ratios 1.25 and 4 separate mild and strong heterogeneity using the same random maps. They do not represent a fitted material. The 2% tolerance describes one user's numerical accuracy budget. The separate decision threshold is the uniform energy, obtained from the reference solve. The 10⁻¹⁰ relative residual is a high-accuracy comparator. Seed blocks 81000–81047 and 91000–92023 separate method timing and frequency estimation. The Monte Carlo standard error pertains to draws under the specified field law.

Numerical interval containment is checked against direct solution to 10⁻⁹ absolute tolerance. Formal floating-point interval arithmetic was not used. The scripts use a safety margin of 10⁻¹¹ when deciding whether to escalate threshold queries. An unresolved equality raises an error rather than silently classifying it. No equality occurred. Machine-specific timings should be regenerated on the hardware of interest. No single-host timing ranks all numerical algorithms for a family of materials.
