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

## Two-size benchmark: data and estimands

* `setup.csv` contains geometry setup shared by all approaches, uniform factorization and right-hand-side solve required by all recycled approaches, extra trial-field setup and extra 16-column Woodbury setup. Direct SuperLU requires geometry only.
* `queries_64.csv` and `queries_208.csv` contain 48 fresh fields per family and size. Solver timings include field-dependent matrix assembly for those that need it. Field generation is outside each compared timer. The fixed-support shortcuts exploit the 16 prescribed edge indices. The dense-rank Woodbury storage estimate is the size of one dense edge-by-edge matrix, not a lower bound on every possible implementation. Direct LU solves the sparse system from scratch each time. The 2% PCG implementation and tight PCG use the same preconditioner and initial guess. The screened hybrid's point-query column sums individual per-field component times, including fallback when required. Independent runs should reproduce numerical conclusions within double-precision tolerance. Latencies and bit patterns can differ.
* `uq_draws.csv` gives the certified interval, decision and optional direct-solver check for every field in the many-query experiment. `uq_summary.csv` separates per-family input generation, actual query pipeline and sentinel-check times from shared amortized setup. Use the threshold decision directly when estimating an exceedance probability. A certified 2% point estimate alone can give the wrong indicator.
* `CITATION-CHECKS.md` identifies the source sections underlying the paper. Timing tests compare sparse LU, known-support Woodbury, uniform-factor PCG and two forms of stopping. Multigrid and continuum approximations have no timing entries in these experiments.

## Stipulated design choices

An even square side length of 64 provides a smaller comparison; 208 provides 42,848 free nodes and 86,059 surviving bonds. Bond spacing, intact tip stiffness and baseline stiffness are one. Both grip rows have prescribed displacement. Sixteen nearest surviving bonds define the known-support case, keeping the modification low rank without assuming all perturbations are weak. Their stiffnesses vary independently and log-uniformly between a quarter and four times unit stiffness. Distributed fields threshold Gaussian noise smoothed over three lattice spacings, producing connected clusters. Ratios 1.25 and 4 separate mild and strong heterogeneity using the same random maps. They do not represent a fitted material. The 2% tolerance describes one user's numerical accuracy budget. The separate decision threshold is the uniform energy, obtained from the reference solve. The 10⁻¹⁰ relative residual is a high-accuracy comparator. Seed blocks 81000–81047 and 91000–92023 separate method timing and frequency estimation. The Monte Carlo standard error pertains to draws under the specified field law.

Numerical interval containment is checked against direct solution to 10⁻⁹ absolute tolerance. Formal floating-point interval arithmetic was not used. The scripts use a safety margin of 10⁻¹¹ when deciding whether to escalate threshold queries. An unresolved equality raises an error rather than silently classifying it. No equality occurred. Machine-specific timings should be regenerated on the hardware of interest. No single-host timing ranks all numerical algorithms for a family of materials.

## Solver and size comparison

The five-size comparison uses the same network and random-field generator in `benchmark.py`. On Linux x86-64 with Python 3.12, create a fresh environment and run from this directory:

```sh
python -m pip install --require-hashes -r scaling-requirements.lock
python scaling.py --out reproduced-scaling
python verify_scaling.py
```

The verifier reads the distributed tables. To inspect timings regenerated on another machine, copy the verifier into the output directory and invoke it there. Timing results depend on machine load. The seeds, matrix properties and responses are portable subject to numerical library variation. `plot_scaling.py` recreates `scaling.png` from `scaling_summary.csv` with Matplotlib; the figure was rendered with the version given below.

The default run uses five even square sizes divisible by four (64, 96, 144, 208, 288) and seed indices 183000–183007. The same eight indices pair the field families, and every method at a given family/size/seed receives exactly the same field. Field correlation scale is three lattice spacings regardless of size. Family and method order are randomized separately within size using seed 247135 + n. A 24-node-side warm-up runs before all reported timings. The results were obtained on Linux x86-64, Intel Xeon Platinum 8573C, in a shared execution environment. Four thread-limit variables are set to one before imports. Fresh numerical assemblies reuse the index structure in `Network`. Values, solver conversion, factorization or hierarchy construction and solution are charged per field. Residual checks, matching-solution checks and memory-array extraction are outside the timers. The two AMG PCGs start from zero. The output pipeline uses the reference solution, which is included in setup.

- `scaling_fields.csv` has six method measurements per field, for 720 records. Fresh SuperLU defaults to SciPy's ordering. Fresh CHOLMOD converts the lower triangle, performs symbolic analysis and numeric factorization and solves. Reused-pattern CHOLMOD only refills its existing sparse matrix, refactorizes numerically and solves. CVXOPT 1.3.3's wheel loads libcholmod.so.3.0.6 and uses supernodal LLᵀ. Fresh and uniform AMG are classical Ruge–Stüben hierarchies, with default strength threshold 0.25, classical interpolation, symmetric pre/post Gauss–Seidel sweeps, a V-cycle and a pseudoinverse on at most 10 coarse variables. The uniform hierarchy's initial V-cycle charges lazy coarse setup to its recorded preprocessing. SciPy CG uses zero absolute tolerance and 1e-10 relative tolerance.
- `output_hybrid` times the actual scalar screen followed, where needed, by PCG preconditioned by a uniform CHOLMOD factor. Its total timer includes internal sparse assembly only when the screen fails. For this method `solve_ms` is the entire pipeline. Zero in its separate setup/assembly columns means included in that pipeline timer, not that fallback assembly is costless. It returns a scalar and therefore has no full-vector residual record.
- `scaling_summary.csv` gives within-cell medians by phase and total, plus min/max totals. A median across sums generally differs from a sum of component medians. `scaling_setup.csv` separates shared geometry and uniform assembly from explicit reusable setup. For cold-start costs, add shared geometry to all methods and, for each reused method, add uniform assembly and its own one-time setup. A single-geometry timing cannot predict all applications.
- `scaling_memory.csv` counts exported numeric/index sparse arrays for seed 183000 per size and family. SuperLU's L and U include diagonals, CVXOPT's extracted L uses its own index type, and AMG counts distinct backing arrays for A, P and R across levels. The last measure includes the fine A. These partial storage measures exclude peak resident memory, workspace, symbolic objects and dense coarse solves.
- `scaling_environment.json` pins the measured solver stacks. No third-party library binaries are distributed; the lockfile points to the upstream wheels.

Figure rendering: Matplotlib 3.10.8.

## Independent energy check

Run `python independent_energy.py` to regenerate `independent_energy.csv`. A separate implementation loops through horizontal and vertical bonds, assembles the full dense stiffness operator and imposes displacement zero on the bottom grip and one on the top. It solves pre-cut and post-cut equilibrium separately and sums individual spring energies. Dividing their difference by the squared tip force gives a direct check on the assembly, boundary conditions and dipole normalization in the production code. The nine cases combine widths 8, 12 and 16 with uniform, smoothly varying and stepped stiffnesses. This check tests the numerical observable within the specified mechanical model.


## Timing spreads, setup recovery and threshold classification

Run the following from this folder. The replay uses the pinned Python 3.12 Linux solver stack in `scaling-requirements.lock`.

```sh
python timing_spread.py
python threshold_compare.py --out reproduced-threshold
python verify_threshold_compare.py
```

`timing_spread.py` reads the distributed 720-row point-solver timings. For each size, family and method, `scaling_quartiles.csv` records the linearly interpolated 25th percentile, median, 75th percentile and interquartile width over the eight simulated stiffness fields. Each received one timed execution per method. Consequently field variation and machine timing noise both enter these descriptive spreads. `scaling_amortization.csv` subtracts direct-symbolic setup from the certificate setup; both require the same uniform assembly and geometry. It gives the number of queries needed to offset the premium assuming constant median charges and the first observed crossing when the eight paired savings are accumulated by increasing seed. A blank field means that no crossing was observed in the specified eight queries, or that the median savings were nonpositive, depending on the column.

`threshold_compare.py` reconstructs all 3,072 fields and checks them against `uq_draws.csv` while it executes. Both uniform factors pay for their preconditioner, trial weights and 16-bond Woodbury cache. The field Cholesky baseline pays once for symbolic analysis, then assembles and numerically refactorizes the field matrix on every query. All three operate on exactly the same stiffness values. Seed 712308 permutes the family–field cases and the three timings within each case; a disposable 32-wide network warms up the paths. The 96 archived SuperLU spot checks give a separate reference for the new direct Cholesky values. The new hardware and versions are recorded in `threshold_environment.json`; compare each method with contemporaneous measurements on the same host.

The distributed `threshold_fields.csv`, `threshold_setup.csv` and `threshold_summary.csv` belong to this paired replay. Timing totals exclude generation, common setup and the comparisons performed outside the timers. `verify_threshold_compare.py` independently recomputes each decision and the summary totals from these files and reads no solver code. To check a regenerated timing run, copy this verifier and `uq_draws.csv` into its output directory and execute the copy. The older `uq_*.csv` files and legacy solver runs are retained for numerical comparison.
