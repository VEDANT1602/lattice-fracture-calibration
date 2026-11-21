# When a stiffness query can stop before a point solve

Vedant Thaker  
Independent Researcher, Jeffersonville, Indiana, USA  
vedant.thaker1602@gmail.com

Numerical methods report

### Abstract

Estimating a response to a specified accuracy and deciding whether it crosses a threshold can require different work. On a cut spring network with 42,848 free nodes, trial energy bounds settled 854 of 1,024 threshold queries for strongly heterogeneous fields. The same bounds met a 2% point-value tolerance on none of 48 strong fields in a separate panel. Sparse Cholesky with a reusable symbolic factor was a substantially faster direct comparator than fresh SuperLU. On eight other strong fields of this size, median query times were 30.7 ms for reused-pattern Cholesky, 138.8 ms for SuperLU and 20.7 ms for a 2% output certificate. Rebuilding a classical algebraic multigrid hierarchy and converging preconditioned conjugate gradients took 161.5 ms. Five specimen sizes establish the range of this comparison. The practical choice depends on the requested output, the contrast and the number of queries that can share setup.

## A value or a threshold decision

Single-bond changes in elastic networks can substantially alter a global response [1]. Here the quantity of interest is local. For a cut of prescribed length and location, how much elastic energy is released per squared pre-cut force in the next bond? One user may ask for that number to 2% relative accuracy. Another may count the fraction of random stiffness fields for which it exceeds the uniform-network value. A 2%-accurate answer near the reference can fall on its wrong side. A wide interval far from the reference can settle the decision immediately.

Consider an n × n square network of scalar, nearest-neighbor springs at unit spacing. The two grip rows have prescribed displacement, and the lateral edges are free. Across the horizontal midline, remove the n/4 leftmost vertical bonds (n divisible by four), then release one further tip bond of unit stiffness. Let K be the positive-definite stiffness matrix of the remaining free nodes after this cut, and h its signed unit dipole at the tip endpoints. Under fixed grips, removing the tip releases its stored energy plus the relaxation energy of the surviving network. For nonzero pre-cut tip force F, these are F²/2 and F²hᵀK⁻¹h/2. Thus the normalized release E equals (1 + R)/2, where R = hᵀK⁻¹h. The unit stiffness fixes the scale. The index n measures specimen size at fixed spacing. A separate dense assembly solved displacement-controlled equilibria before and after removing the tip. Across three small sizes and three stiffness patterns, the measured energy loss divided by the squared tip force agreed with the dipole calculation within 10⁻¹¹ relative error.

## What a uniform solve buys

Preprocessing can pay off when thousands of parameter queries share a geometry [2]. Factor the uniform post-cut matrix K₀ once. Write x₀ = K₀⁻¹h and R₀ = hᵀx₀. For bond incidence bₑ, the nonnegative weights wₑ = (bₑᵀx₀)²/R₀ sum to one. Use the uniform-network potential and flow as admissible trial fields in the Dirichlet and Thomson principles [3]. For stiffness vector k this gives (1 + R₀/Σₑwₑkₑ)/2 ≤ E ≤ (1 + R₀Σₑwₑ/kₑ)/2. For 0 < L ≤ E ≤ U, returning √(LU) incurs at most √(U/L) − 1 relative error. To classify E, compare both endpoints with the reference.

When only 16 known bonds change, cache H = BᵀK₀⁻¹B and f = Bᵀx₀ for their incidence columns B. With diagonal stiffness changes Δ, Woodbury gives R = R₀ − fᵀ(I + ΔH)⁻¹Δf [4]. Rank one is the Sherman–Morrison special case. For a distributed change on 86,059 bonds, a dense matrix with one double-precision entry per changed-bond pair alone would occupy about 55 GiB. The benchmark uses this classical update only for the stipulated 16-bond support.

For general fields, start preconditioned conjugate gradient (PCG) at x₀ and apply K₀⁻¹ as preconditioner [5]. One variant targets relative residual 10⁻¹⁰. Another stops for the requested scalar output. For an iterate x set r = h − Kx and L = (1 + 2hᵀx − xᵀKx)/2. The unaccounted energy equals rᵀK⁻¹r/2. Since K ≥ kₘᵢₙ K₀ in quadratic-form order, U = L + rᵀK₀⁻¹r/(2kₘᵢₙ) bounds E from above. The stopping test recomputes the true residual before recording the interval. This certificate assumes exact arithmetic. The numerical checks below use a tolerance of 10⁻⁹.

## A direct comparison on identical fields

The local family varies the 16 nearest surviving bonds to the tip independently, log-uniformly between 1/4 and 4. Distributed fields start with unit normal noise on the nodes, smoothed at three lattice spacings. Each bond takes its high or low stiffness according to the sign of the sum at its endpoints. Paired maps give mild values 0.8 and 1.25, or strong values 0.25 and 4. These laws define computational loads. No material has been fitted to them.

Five sizes have 3,968 to 82,368 free nodes. At each size, eight seed indices generated fields in each family. Six methods ran in randomized order on each field. The two direct CHOLMOD configurations use supernodal LLᵀ [6]. One orders and factorizes each field from scratch. The other retains its symbolic factor and lower-triangle layout, then replaces the values and refactorizes. Fresh SuperLU uses its default sparse ordering. Classical Ruge–Stüben algebraic multigrid (AMG) in PyAMG [7] preconditions CG, either through a newly built hierarchy per field or through one hierarchy built on uniform stiffness. Both use symmetric Gauss–Seidel smoothing, V-cycles and a relative residual tolerance of 10⁻¹⁰. Those five methods solve for the full displacement vector.

The sixth method only guarantees the scalar response to 2%. It tries the trial-field interval before assembly, then, when necessary, advances PCG using a uniform Cholesky factor until the residual interval is narrow enough. Both iterative multigrid methods start from zero; the certified iteration starts from the uniform solution paid for during setup. Single-thread timings include per-field assembly, required representation changes, factorization or hierarchy construction, and solution. Timing excludes field generation, shared geometry and reusable setup. Stored fields and all individual timings accompany the paper.

At n = 208, common geometry cost 19.2 ms. Uniform assembly added 1.6 ms for the reusable methods. Further setup costs were 19.9 ms for the Cholesky pattern, 41.6 ms for uniform AMG, and 52.3 ms for the Cholesky-based certificate. On strong fields, preserving the Cholesky pattern cut the direct median from 50.9 to 30.7 ms. The scalar-only method took 20.7 ms, a factor of 1.48 relative to that direct baseline. Comparing it only to the 138.8 ms SuperLU result would imply a factor of 6.71. The repeated-query benefit is smaller once the direct solver uses symmetry, and the scalar method’s larger setup matters for short batches. At the same size the mild-field certificate took 0.17 ms because the trial interval sufficed before assembly.

| Side n | Free nodes | Fresh LU | Fresh Cholesky | Reused Cholesky | Field AMG | Uniform AMG | 2% scalar |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 64 | 3,968 | 8.1 | 4.0 | 1.8 | 14.9 | 29.3 | 1.6 |
| 96 | 9,024 | 18.6 | 9.2 | 4.5 | 31.9 | 57.9 | 3.9 |
| 144 | 20,448 | 51.7 | 20.6 | 12.0 | 69.8 | 115.8 | 10.4 |
| 208 | 42,848 | 138.8 | 50.9 | 30.7 | 161.5 | 228.4 | 20.7 |
| 288 | 82,368 | 337.8 | 110.3 | 69.7 | 329.4 | 424.2 | 44.4 |
Table 1. Median milliseconds per query across eight strong fields per size. Field AMG builds the hierarchy anew. Uniform AMG reuses its hierarchy. The scalar method has a different accuracy requirement.

![Scaling of median per-field wall time for the six methods](study/scaling.png)

Figure 1. Strong distributed fields across five specimen sizes. Time is on a logarithmic scale. Each solid curve measures a full solve. The dashed curve measures the time to a certified scalar estimate.

Of the 40 strong-field comparisons in Table 1, the certified result was faster than reused-pattern Cholesky in 35, including all eight at each of the two largest sizes. The iterative AMG methods cost more than reused-pattern Cholesky at every size on median. At n = 288, building a fresh hierarchy reduced median CG iterations from 46.5 to 24 relative to uniform reuse. The medians of the phase times and individual paired timings are available in the supplement. On all 120 fields, the five full-solve responses differed from SuperLU by at most 1.9 × 10⁻¹⁴ in relative terms. The largest observed certified-estimate error was 1.478%.

The storage tradeoff differs from the latency ranking. For one strong field at n = 288 the exported Cholesky L arrays occupied 53.6 MiB, the exported SuperLU L and U arrays 90.0 MiB, and the sparse operators and transfers in the nine AMG levels 18.4 MiB. These array totals use the actual index widths of each binding. Process workspaces, symbolic data, transient allocations and coarse dense factors are excluded. They do not establish peak memory use.

## Classifying a threshold

A 48-field-per-family accuracy panel at each of n = 64 and 208 found that the trial interval alone met the 2% point-value criterion on all mild distributed fields and on no local or strong fields. On the local fields, the known-support Woodbury update matched the direct value to roundoff. This rank-16 calculation illustrates the benefit of knowing which bonds vary; the update identity is classical [4].

At n = 208 the uniform threshold is E₀ = 1.2059805. For 1,024 additional seed indices in each family, the localized queries used Woodbury. Distributed queries invoked the certified iteration only if the uniform trial bounds touched the threshold. On 32 seed indices per family, direct solves independently checked the classifications; all 96 agreed.

| Field family | Above E₀ / 1,024 | Fraction | Iterative fallbacks | Query time (s) |
| --- | --- | --- | --- | --- |
| 16 local | 555 | 0.542 | 0 | 0.042 |
| Distributed mild | 508 | 0.496 | 33 | 0.931 |
| Distributed strong | 532 | 0.520 | 170 | 4.857 |
Table 2. Timed threshold pipeline over 1,024 fields per family. The stored threshold implementation uses a reusable SuperLU factor for the uniform system. Query time excludes generation and 0.349 s shared setup. This is a separate timing block from Table 1.

On strong fields the trial interval settled 854 threshold decisions without iteration. In contrast, the same criterion settled none of the strong point estimates in the 48-field accuracy panel. Width controls the value request; distance from the reference controls the threshold request. The Monte Carlo standard error for the modeled exceedance fraction is about 0.016. Acceptance fractions depend on the field law, cut geometry and reference.

## Limits of the timing comparison

These are finite, two-dimensional scalar networks with fixed physical spacing, a prescribed cut and one load vector. Timing on a shared host, library versions, ordering, implementation overhead and multigrid tuning all affect the measurements. Five sizes show the observed scaling range. They cannot establish an asymptotic complexity law, a faster solver on three-dimensional systems, or a prediction of physical crack growth. A user who needs the full solution cannot substitute the scalar certificate.

## Data and code availability

The supplement includes scripts, environment locks, source checks and machine-readable results. The five-size comparison has 720 timed method–field records. The separate accuracy and threshold datasets contain 288 and 3,072 field records. All fields can be regenerated from the published seeds.

[Study files, software and source records](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/main/study)

## Author and disclosures

Vedant Thaker is the sole author, set the solver-selection question and requested direct and multigrid comparators at larger scales. OpenAI Codex assisted with the derivations, implementation, benchmark execution, numerical cross-checks and manuscript drafting.

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors. The author declares no competing interests.

## References

[1] D. Hexner, A. J. Liu and S. R. Nagel. Linking microscopic and macroscopic response in disordered solids. Physical Review E 97, 063001 (2018). <https://doi.org/10.1103/PhysRevE.97.063001>

[2] G. Rozza, D. B. P. Huynh and A. T. Patera. Reduced Basis Approximation and a Posteriori Error Estimation for Affinely Parametrized Elliptic Coercive Partial Differential Equations: Application to Transport and Continuum Mechanics. Archives of Computational Methods in Engineering 15, 229–275 (2008). <https://doi.org/10.1007/s11831-008-9019-9>

[3] A. Nachmias. Random Walks and Electric Networks. In Planar Maps, Random Walks and Circle Packing, Lecture Notes in Mathematics 2243, 11–31 (2020). <https://doi.org/10.1007/978-3-030-27968-4_2>

[4] H. V. Henderson and S. R. Searle. On Deriving the Inverse of a Sum of Matrices. SIAM Review 23(1), 53–60 (1981). <https://doi.org/10.1137/1023004>

[5] R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra, V. Eijkhout, R. Pozo, C. Romine and H. van der Vorst. Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods. SIAM (1994). <https://doi.org/10.1137/1.9781611971538>

[6] Y. Chen, T. A. Davis, W. W. Hager and S. Rajamanickam. Algorithm 887: CHOLMOD, Supernodal Sparse Cholesky Factorization and Update/Downdate. ACM Transactions on Mathematical Software 35(3), article 22 (2008). <https://doi.org/10.1145/1391989.1391995>

[7] N. Bell, L. N. Olson and J. Schroder. PyAMG: Algebraic Multigrid Solvers in Python. Journal of Open Source Software 7(72), 4142 (2022). <https://doi.org/10.21105/joss.04142>
