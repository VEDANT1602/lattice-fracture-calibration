# Choosing a solver for repeated stiffness queries on a cut spring network

Vedant Thaker  
Independent Researcher, Jeffersonville, Indiana, USA  
vedant.thaker1602@gmail.com

21 November 2026

### Abstract

A stiffness-field study may require the value of a local response or just the fraction of fields exceeding a reference. Those requests need different stopping rules. A scalar square network with 42,848 unknowns provides a repeated-query benchmark. When 16 specified springs vary, a reusable Woodbury update gives the direct-solve result in 0.14 ms per field after setup. A trial-field bound resolves a 2% relative-accuracy request for all 48 mildly heterogeneous distributed fields tested and none of 48 strongly heterogeneous fields. Even at strong contrast, the same bound classifies 854 of 1,024 fresh fields relative to a uniform reference without iteration. The location of the changing springs and the requested output both affect solver choice.

## A value or a threshold decision

Single-bond changes in elastic networks can substantially alter a global response [1]. Here the quantity of interest is local. For a cut of prescribed length and location, how much elastic energy is released per squared pre-cut force in the next bond? One user may ask for that number to 2% relative accuracy. Another may count the fraction of random stiffness fields for which it exceeds the uniform-network value. A 2%-accurate answer near the reference can fall on its wrong side. A wide interval far from the reference can settle the decision immediately.

Consider an n × n square network of scalar, nearest-neighbor springs at unit spacing. The two grip rows have prescribed displacement, and the lateral edges are free. Across the horizontal midline, remove the n/4 leftmost vertical bonds, then release one further tip bond of unit stiffness. Let K be the positive-definite stiffness matrix of the remaining free nodes after this cut, and h its signed unit dipole at the tip endpoints. Under fixed grips, removing the tip releases its stored energy plus the relaxation energy of the surviving network. For nonzero pre-cut tip force F, these are F²/2 and F²hᵀK⁻¹h/2. Thus the normalized release E equals (1 + R)/2, where R = hᵀK⁻¹h. The unit stiffness fixes the scale. The index n measures specimen size at fixed spacing.

## What a uniform solve buys

Preprocessing can pay off when thousands of parameter queries share a geometry [2]. Factor the uniform post-cut matrix K₀ once. Write x₀ = K₀⁻¹h and R₀ = hᵀx₀. For bond incidence bₑ, the nonnegative weights wₑ = (bₑᵀx₀)²/R₀ sum to one. Use the uniform-network potential and flow as admissible trial fields in the Dirichlet and Thomson principles [3]. For stiffness vector k this gives (1 + R₀/Σₑwₑkₑ)/2 ≤ E ≤ (1 + R₀Σₑwₑ/kₑ)/2. For 0 < L ≤ E ≤ U, returning √(LU) incurs at most √(U/L) − 1 relative error. To classify E, compare both endpoints with the reference.

When only 16 known bonds change, cache H = BᵀK₀⁻¹B and f = Bᵀx₀ for their incidence columns B. With diagonal stiffness changes Δ, Woodbury gives R = R₀ − fᵀ(I + ΔH)⁻¹Δf [4]. Rank one is the Sherman–Morrison special case. For a distributed change on 86,059 bonds, a dense matrix with one double-precision entry per changed-bond pair alone would occupy about 55 GiB. We limit this implementation to the localized workload.

For general fields, start preconditioned conjugate gradient (PCG) at x₀ and apply K₀⁻¹ as preconditioner [5]. One variant targets relative residual 10⁻¹⁰. Another stops for the requested scalar output. For an iterate x set r = h − Kx and L = (1 + 2hᵀx − xᵀKx)/2. The unaccounted energy equals rᵀK⁻¹r/2. Since K ≥ kₘᵢₙ K₀ in quadratic-form order, U = L + rᵀK₀⁻¹r/(2kₘᵢₙ) bounds E from above. We recompute the true residual before recording the interval. This certificate assumes exact arithmetic. The numerical checks below use a tolerance of 10⁻⁹.

## Held-out timing comparison

The benchmark uses 48 independently seeded fields per family and size. In the local family, stiffnesses on the 16 bonds nearest the tip are independent and log-uniform from 1/4 to 4. Distributed fields threshold nodewise Gaussian white noise after Gaussian smoothing with scale three lattice spacings, using the sign of the sum at each bond’s endpoints. The two contrasts share field patterns. Mild stiffnesses take values 0.8 or 1.25, and strong ones take 0.25 or 4. These choices distinguish small and large departures without fitting any material. Preprocessing uses geometry, stipulated support and the uniform network alone.

Runs used one thread, randomized query and method order, precomputed sparsity, and a fresh SuperLU factorization as reference on every field. Timings include field-dependent numeric assembly when needed. Field generation is excluded equally. At n = 208, geometry cost 38.1 ms for all methods. The reusable factor and reference solve cost 161.6 ms. The bound required a further 0.75 ms to form weights, and Woodbury required 40.5 ms for its 16 response columns and small matrix. Table 1 gives means across fields on this host. Screen-plus-PCG sums separately timed components for each case. The next experiment times a directly executed pipeline.

The screen resolved 0, 48 and 0 of the 48 point-value requests in the local, distributed mild and distributed strong families. The same counts occurred at n = 64 (3,968 unknowns). At strong contrast, fully converging PCG cost more than refactorization at both sizes. Stopping PCG as soon as it certified the output made it faster here. On the local fields, Woodbury recovered the direct values to floating-point accuracy. One local field at n = 208 changed threshold classification when its permitted 2% point estimate replaced its direct value.


| Field family | Direct LU | Tight PCG | 2% PCG | Screen + 2% PCG | Woodbury |
| --- | --- | --- | --- | --- | --- |
| 16 local | 185.4 | 130.2 | 34.32 | 34.35 | 0.142 |
| Distributed mild | 183.1 | 99.4 | 21.86 | 0.160 | — |
| Distributed strong | 199.5 | 390.2 | 44.93 | 45.08 | — |
Table 1. Mean milliseconds per field at n = 208, with 48 held-out fields per row. Setup is reported separately. A dash marks a method omitted at distributed rank.

## Exceedance frequency under three input laws

At n = 208, set the threshold to the uniform response, E₀ = 1.2059805. Draw 1,024 fresh fields per family. Use exact Woodbury for known support. For distributed fields, apply the uniform-field bounds and switch to certified PCG only when the interval touches the reference. The timings in Table 2 include this control flow. On 32 seed indices selected in advance per family, a fresh direct solve checked the classifications. All 96 agreed.

At strong contrast the bounds settled 854 of 1,024 decisions on their own, despite settling no point-value requests at the 2% tolerance in Table 1. The estimated exceedance fractions have Monte Carlo standard error near 0.016 under each chosen field law. These are separate model probabilities. The two distributed rows use corresponding random maps, so their estimates are paired. All queries together reused 0.349 s of setup in this second run. The pattern generalizes algebraically to other positive networks, but the acceptance rates and timings depend on the field law, topology, solver and hardware. The benchmark measures a prescribed scalar cut on two finite grids. Extrapolation to crack growth, material failure or measured toughness requires different evidence.


| Field family | Above E₀ / 1,024 | Fraction | PCG fallbacks | Generate (s) | Queries (s) |
| --- | --- | --- | --- | --- | --- |
| 16 local | 555 | 0.542 | 0 | 0.080 | 0.042 |
| Distributed mild | 508 | 0.496 | 33 | 3.172 | 0.931 |
| Distributed strong | 532 | 0.520 | 170 | 3.157 | 4.857 |
Table 2. Each row uses 1,024 fields. Cumulative wall times exclude the independent direct-solver checks. Fixed-support queries use Woodbury throughout.

## Data and code availability

The companion supplement includes benchmark.py, uncertainty.py, the environment lock, 288 field-level timing records, 3,072 threshold decisions and the verification script. Each random field is regenerated from its recorded seed.

[Study files, software and source records](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/main/study)

## Author and disclosures

Vedant Thaker is the sole author and set the research question. OpenAI Codex assisted with deriving and implementing the numerical comparisons, executing the random-field benchmarks and drafting the results.

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors. The author declares no competing interests.

## References

[1] D. Hexner, A. J. Liu and S. R. Nagel. Linking microscopic and macroscopic response in disordered solids. Physical Review E 97, 063001 (2018). <https://doi.org/10.1103/PhysRevE.97.063001>

[2] G. Rozza, D. B. P. Huynh and A. T. Patera. Reduced Basis Approximation and a Posteriori Error Estimation for Affinely Parametrized Elliptic Coercive Partial Differential Equations: Application to Transport and Continuum Mechanics. Archives of Computational Methods in Engineering 15, 229–275 (2008). <https://doi.org/10.1007/s11831-008-9019-9>

[3] A. Nachmias. Random Walks and Electric Networks. In Planar Maps, Random Walks and Circle Packing, Lecture Notes in Mathematics 2243, 11–31 (2020). <https://doi.org/10.1007/978-3-030-27968-4_2>

[4] H. V. Henderson and S. R. Searle. On Deriving the Inverse of a Sum of Matrices. SIAM Review 23(1), 53–60 (1981). <https://doi.org/10.1137/1023004>

[5] R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout, R. Pozo, C. Romine and H. van der Vorst. Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods. SIAM (1994). <https://doi.org/10.1137/1.9781611971538>
