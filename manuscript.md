# When one reference solve can replace a new solve for a cut spring network

**Vedant Thaker**  
Independent Researcher, Jeffersonville, Indiana, United States  
Correspondence: vthaker1@asu.edu  
21 November 2026

## Abstract

An isolated query for the energy released by a prescribed spring cut needs only one sparse solve. A reduced calculation has to recover its setup cost when many stiffness fields share the same lattice and cut. One uniform-network solve supplies trial forces and displacements that bracket the release for every positive stiffness field on that topology. On a panel of 2,160 fields covering three patch shapes and two random-field models, a 2% relative-error threshold accepted 1,037 predictions without a new solve. The remaining queries used the exact solver. The largest observed error of the combined procedure was 1.285%. In mixed batches the combined calculation was 1.79–2.03 times as fast after calibration as solving every field. None of the 288 strongest random-field cases passed the bound. Direct solution remains the appropriate choice for that batch and for an isolated query. The experiment concerns energy after a specified cut. It does not predict which bond will break.

## 1. When is a calibration worth keeping?

Removing a loaded bond changes the displacements available to the rest of a network. Hexner, Liu and Nagel related local bond response to changes in global elastic moduli after bond removal [1]. Reduced-basis methods have long separated offline calibration from repeated evaluation and used output bounds to assess accuracy [2]. Here the output is a single energy difference for a fixed cut. The question is how often one homogeneous solution can supply an acceptable estimate when the location, shape and values of the stiffness changes are unknown at calibration.

An earlier circular-patch calculation compared 288 nonuniform configurations. A geometry-matched uniform multiplier gave 4.18% mean absolute relative error, and a local slope correction gave 0.43%. That slope needed two perturbed solves for each patch description. Fitting a separate slope to each new patch is a poor bargain if a single direct solve already answers the question. The present test sets the shape-specific method aside and calibrates once per specimen and cut. The unequal costs of calibration, screening and fallback are measured separately.

## 2. The exact one-solve baseline

Each spring connects horizontal or vertical neighbours on an n × n square lattice with one scalar displacement per node and unit spacing. The top and bottom rows are grips. The side boundaries are free. Starting at the left edge, remove a vertical springs that cross the horizontal midline. The next vertical spring is the prescribed cut. At a unit force magnitude in that spring, its stored energy is 1/(2kₜ), where kₜ is its stiffness. The grip positions stay fixed during subsequent relaxation.

Let K be the stiffness matrix of the network after the specified bond is cut, restricted to free displacements. Let h apply a unit positive force to one former endpoint and a unit negative force to the other. Releasing the bond unloads this dipole. The subsequent relaxation contributes half the compliance R = hᵀK⁻¹h. The exact release is D = (1/kₜ + R)/2. Positive stiffness and connection to the grips make this solve well posed for all tested fields. The direct comparator assembles the cut matrix, factors it once and solves one right-hand side. It requires neither an intact-network solve nor subtraction of two large elastic energies.

The intact and cut calculations also agree by a classical rank-one identity. Sherman and Morrison considered a change to a single matrix element in 1950 [3]. Henderson and Searle give the general outer-product update and its historical development [4]. With intact matrix A = K + kₜhhᵀ, let c = hᵀA⁻¹h. Applying the Sherman–Morrison formula to K = A − kₜhhᵀ gives R = c/(1 − kₜc), hence D = 1/[2kₜ(1 − kₜc)]. This establishes equivalence with the intact-network expression. The timed baseline solves on K directly.

## 3. A bracket from the uniform network

Solve K₀u₀ = h once with unit background stiffness, keeping the topology and grips fixed. Set R₀ = hᵀu₀. For each remaining bond e, let qₑ be its squared displacement difference under u₀, divided by R₀. These nonnegative weights sum to one. Record zero weight for the cut bond. For a new field, all bond stiffnesses kₑ and the separate tip stiffness kₜ must be supplied.

The uniform dipole flow remains a feasible unit flow for any positive kₑ. Its complementary energy bounds R above by R₀ times the weighted arithmetic mean of 1/kₑ. Applying the unit-voltage trial u₀/R₀ to the new network bounds R below by R₀ divided by the weighted arithmetic mean of kₑ. These are the finite-network Thomson and Dirichlet variational principles [5]. The fixed grips become one zero-increment-displacement boundary node. The bounds therefore apply to the scalar incremental problem even though the initial intact loading had nonzero grip separation.

Write L and U for the release limits obtained by substituting the two compliance limits into D. The screen returns C = √(LU) if √(U/L) − 1 ≤ ε. Otherwise it performs the exact solve. For D anywhere in the interval, the maximum absolute relative error of C cannot exceed √(U/L) − 1. This certificate follows from the energy inequalities. It does not use labels from the test panel. The estimate costs two weighted sums and one reciprocal per edge with no factorization.

Two point estimates supply controls: keeping the alternative-path compliance at R₀, and approximating it by R₀ exp(−Σ qₑ log kₑ). The latter matches the first variation of log R at unit stiffness: differentiating hᵀK⁻¹h gives a negative squared bond displacement for each stiffness derivative. Both point estimates retain the exact tip term. Their errors show how much the fallback contributes.

## 4. Test fields and error denominators

The six fixed topologies cross n = 16, 32, 48 with a = n/4 or n/2. Each has a distinct uniform calibration and 360 test queries. No nonuniform field is used to fit the weights. Three compact masks are a rotated ellipse (semiaxes 3.1 and 1.45), a rotated strip (half-length 5, half-width 0.85), and two disjoint disks (radii 1.4 and 1.8, centers separated by 4.8). Each mask is translated within three lattice spacings of the tip and rotated through a seeded angle. These geometries test a broad inclusion, a slender inclusion, and disconnected support. They are prescribed stress tests, with dimensions in lattice spacings.

For local masks, covered bonds have stiffness r or 1/r and background bonds have unit stiffness. For distributed fields, every bond has either r or 1/r. The independent binary model uses fair independent signs. The correlated binary model smooths nodewise Gaussian noise with a two-spacing Gaussian kernel, bilinearly samples it at bond midpoints and thresholds at the median. Sign reversal exchanges stiff and soft material. The absolute levels are the same across the two binary models. At r = 4 their maximum stiffness ratio is 16, compared with 4 for a local patch.

Each family uses 12 seeds, r = 1.25, 2, 4 and both signs. The contrast settings span mild, intermediate and wide departures from the unit reference. Sharing seeds across geometries and signs provides paired comparisons. The 2,160 configurations are an equally weighted design panel with no distributional claim for physical materials. An error percentage is 100 times |estimate/exact − 1|. A solve is avoided only when the two-sided certificate passes.

## 5. Accuracy and avoided solves

Without fallback the uniform, log-sensitivity and bracket-center estimates had mean errors of 16.03%, 2.22% and 3.22%. The log-sensitivity point estimate had a maximum error of 47.37%, and the bracket center had a maximum error of 57.03%. The worst center prediction was 0.61277 against 0.39023 exactly in a sign-reversed correlated field at n = 16, a = 4 and r = 4 (seed 7). Its interval, 0.35903 to 1.04585, immediately disqualifies the estimate at 2%.

The 2% screen accepted every r = 1.25 case (720), 265 of 720 at r = 2 and 52 of 720 at r = 4. Both binary models required fallback on all 288 strong-contrast cases. Across the full panel, 1,037 queries passed, and the combined method had mean error 0.094% and maximum error 1.285%. At tolerance 0.5%, 1%, 5% and 10%, it avoided 564, 772, 1,443 and 1,700 solves, with maximum observed errors of 0.336%, 0.504%, 3.519% and 8.048%, respectively. The 2% choice is an illustrative numerical specification, not a measured industrial tolerance.

| Field family | Cases | Accepted at 2% | Mean raw error (%) | Max raw error (%) |
| --- | --- | --- | --- | --- |
| Ellipse | 432 | 222 | 1.46 | 14.38 |
| Strip | 432 | 270 | 1.21 | 11.38 |
| Two disks | 432 | 248 | 1.59 | 14.51 |
| Independent binary | 432 | 144 | 5.54 | 54.23 |
| Correlated binary | 432 | 153 | 6.31 | 57.03 |

Table 1. Raw errors are for the geometric center before fallback. Each numerator counts estimates accepted without solving a new system.

## 6. Cost of screening and calibration

The timing test used precomputed stiffness arrays, identical input order within each pair of methods and the same precomputed sparse assembly pattern. A fixed random sample of 90 of the 360 queries for each topology was permuted within each of five repeats. The execution order of the direct, screen-only and hybrid methods was randomized. Reported times are median batch times divided by 90. Exact fallback used the same sparse assembly and SuperLU solver as the timed direct method. Allocation and factorization are included. Shared topology construction, field generation and disk I/O are excluded. One CPU thread was used for numerical libraries on an AMD EPYC 7763.

The screen alone took 0.0067–0.0122 ms per query. With fallbacks included, the speedup over direct solution was 1.79–2.03 times on the six mixed batches. The warm calibration cost is the median of nine new uniform solves with weight construction. The first observed setup time at each topology gives a more conservative amortization check: using those measured costs gives 3–29 queries across topologies, including a single unusually expensive initial setup. The counts assume that the future batch has the same average time savings. Actual break-even for a particular ordering can differ.

A second timing panel consisted entirely of r = 4 independent and correlated binary fields. All 48 queries per topology required fallback. Across topologies, the median hybrid/direct time ratio ranged from 0.94 to 1.11; fluctuations around one offer no consistent advantage. For a batch with that composition, skip calibration and use direct solves.

| n / a | Setup (ms) | Direct (ms) | Hybrid (ms) | Queries to repay setup |
| --- | --- | --- | --- | --- |
| 16 / 4 | 0.603 | 0.400 | 0.223 | 4 |
| 16 / 8 | 0.559 | 0.312 | 0.154 | 4 |
| 32 / 8 | 1.893 | 1.779 | 0.966 | 3 |
| 32 / 16 | 1.798 | 1.509 | 0.789 | 3 |
| 48 / 12 | 4.379 | 4.266 | 2.312 | 3 |
| 48 / 24 | 4.310 | 3.745 | 1.879 | 3 |

Table 2. Time per query for the mixed panel. Break-even rounds setup/(direct − hybrid) up to a whole query. These are measured implementation costs on one machine.

## 7. Scope of the result

The positive result is limited to repeated queries on a known topology, with available stiffness fields, an admissible numerical error and a mixture that includes many fields close enough to the uniform state. A change of cut location, boundary conditions or specimen topology requires a new calibration. The unguarded correction fails as a general replacement for solving the network. Neither the screen nor direct comparison here competes against the best possible reused symbolic factorization, preconditioned iterative solver or higher-dimensional reduced basis.

The energy after a specified cut also differs from the work or path of physical propagation. In the earlier circular grid, a competing spring would reach the same unit threshold first in 23 of 288 nonuniform settings. Hartquist and colleagues combined simulations and experiments to investigate intrinsic energy for crack growth in stretchable networks [6]. Nieves et al. varied directional stiffness and crack speed to calculate a dynamic energy-release-rate ratio for a propagating Mode I crack in a triangular lattice [7]. The present scalar model specifies exactly which bond is cut and assigns the same unit force magnitude at that bond in every trial. It provides a numerical benchmark for energy queries, without a material toughness interpretation.

## Data and code availability

The baseline repository contains the circular-patch reference code and 336-row dataset. The new scripts and machine-readable outputs accompany this paper in the revision supplement. Use transfer.py for the calibration and mixed timings, stress_timing.py for the all-fallback timing experiment, and verify_transfer.py for computational checks. The README and pinned requirements specify the invocation and environment. All fields are generated from recorded seeds.

[Pinned baseline](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/f943322c44df01ecfcb4f81c1b686315a79a9e6a) · [Revision code and data](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/main/revision/supplement) · [Preparation record](https://github.com/VEDANT1602/lattice-fracture-calibration/blob/main/revision/PREPARATION.md)

## Computational verification

All 2,160 releases lay within their independently computed bounds to numerical precision. In 30 sentinel cases the direct dipole release agreed with a separately assembled before/after energy calculation to within 1.27 × 10⁻¹³ relative error. A verifier checks case counts, bound ordering, the five acceptance counts and timing denominators. The previous independent reconstruction of the circular grid also agrees with the pinned release, grip and predictor fields to within 2.6 × 10⁻¹² relative discrepancy across all 336 rows. These checks address the discrete model.

## Author contribution

Vedant Thaker is the sole author. He designed the original scalar lattice fracture study, wrote its manuscript, developed its original software, and requested the comparisons and checks reported in this revision.

## Use of computational assistance

OpenAI Codex assisted with the rank-one and energy-bound calculations, held-out field generation, direct-solve checks, screening and fallback timing comparisons, and reference verification, as well as literature checks and preparation of the revised text. The author set the research question and requested the revision.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

## Competing interests

The author declares no competing interests.

## References

[1] D. Hexner, A. J. Liu and S. R. Nagel. Linking microscopic and macroscopic response in disordered solids. Physical Review E 97, 063001 (2018). https://doi.org/10.1103/PhysRevE.97.063001

[2] G. Rozza, D. B. P. Huynh and A. T. Patera. Reduced Basis Approximation and a Posteriori Error Estimation for Affinely Parametrized Elliptic Coercive Partial Differential Equations: Application to Transport and Continuum Mechanics. Archives of Computational Methods in Engineering 15, 229–275 (2008). https://doi.org/10.1007/s11831-008-9019-9

[3] J. Sherman and W. J. Morrison. Adjustment of an Inverse Matrix Corresponding to a Change in One Element of a Given Matrix. Annals of Mathematical Statistics 21(1), 124–127 (1950). https://doi.org/10.1214/aoms/1177729893

[4] H. V. Henderson and S. R. Searle. On Deriving the Inverse of a Sum of Matrices. SIAM Review 23(1), 53–60 (1981). https://doi.org/10.1137/1023004

[5] A. Nachmias. Random Walks and Electric Networks. In Planar Maps, Random Walks and Circle Packing, Lecture Notes in Mathematics 2243, 11–31 (2020). https://doi.org/10.1007/978-3-030-27968-4_2

[6] C. Hartquist, S. Wang, Q. Cui, W. Matusik, B. Deng and X. Zhao. Scaling Law for Intrinsic Fracture Energy of Diverse Stretchable Networks. Physical Review X 15, 011002 (2025). https://doi.org/10.1103/PhysRevX.15.011002

[7] M. J. Nieves, A. B. Movchan, I. S. Jones and G. S. Mishuris. Propagation of Slepyan’s crack in a non-uniform elastic lattice. Journal of the Mechanics and Physics of Solids 61(6), 1464–1488 (2013). https://doi.org/10.1016/j.jmps.2012.12.006
