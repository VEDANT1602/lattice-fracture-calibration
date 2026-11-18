# How far can a uniform-lattice calibration travel across a stiffness patch?

**Vedant Thaker**  
Working paper · 18 November 2026

## Abstract

Removing one spring from a loaded network releases energy from the surrounding lattice as it relaxes. How well does a calibration from a uniform specimen predict that release near a soft or stiff patch? We prescribe one notch advance in a square scalar-displacement lattice under fixed grips. Across 288 nonuniform configurations, using only the energy stored in the removed spring gives a mean absolute relative error of 57.66%. A uniform calibration matched to the specimen geometry reduces it to 4.18%. A first-order correction for patch stiffness reduces it to 0.43%, with a maximum error of 5.96%. For the tested radii, a small patch centered on the notch tip produces more transfer error than a larger concentric one. Direct energy differences and a rank-one matrix calculation agree to within 6.2e-14 relative error across the full grid. Applying this result to physical fracture energy requires measuring crack propagation.

## One bond, two energies

At a specified breaking force, the energy stored in a linear spring of known stiffness is easy to calculate. Cutting that spring changes which displacements the rest of the network can support. The remaining springs relax, even if the grips never move. A model that assigns the entire energy change to the failed spring omits this relaxation.

Hartquist and colleagues studied homogeneous networks with experiments and simulations across a range of strand behaviours and topologies. Their scaling law relates intrinsic fracture energy per strand areal density to the product of breaking force and stretched strand length. It treats release from the surrounding network as physically consequential [1]. Nieves and colleagues showed that stiffness contrast can change energy release and stability in a moving crack in a triangular lattice [2]. Neither observation prescribes a correction for placing a circular stiffness patch near a finite notch in a clamped scalar lattice.

This computational study takes a narrower measurement than either paper. It holds the notch advance and loading mode fixed, then changes only the stiffness of springs whose midpoints lie inside a patch. Its question is whether knowing the change in local rupture energy suffices to carry a homogeneous calibration to a heterogeneous specimen. A second question is whether two small perturbations near the uniform state can anticipate the sign and magnitude of the error at larger contrast.

## Geometry and reference calculation

The strip is a square of n by n nodes with unit spacing and one scalar displacement per node, interpreted as antiplane motion. Nearest neighbours share a linear spring. The top and bottom rows are held at distinct uniform displacements. Both side boundaries are free. A notch starts at the left edge and removes the first a vertical springs crossing the horizontal midline. The specified advance removes the next vertical spring along that line. This construction keeps the crack increment equal to one bond for all comparisons.

Background stiffness and critical force are both one in dimensionless units. A circular patch multiplies spring stiffness by r. Its midpoint is on the notch midline, with its horizontal coordinate displaced by -4, 0 or +4 lattice spacings from the tip bond. Negative offsets place the patch behind the tip. Two radii, 1.25 and 2.75 spacings, resolve compact and broader inclusions without ambiguous membership at the boundary.

For each case, a sparse equilibrium solve gives the displacement field at unit grip separation. Let f be the magnitude of the resulting force in the tip spring. Linearity gives 1/f as the grip separation needed for the prescribed spring to reach critical force. At this separation the spring contains W = 1/(2k), where k is its stiffness. After cutting the tip spring, the network is equilibrated again with the same grip separation. The total potential energy released, D, is the difference between the two elastic energies. The grips perform no work during the cut. This test specifies which bond breaks and does not check whether some other bond would have failed earlier.

An independent expression checks the direct calculation. Let K denote the intact stiffness matrix restricted to free displacements, and let h map those displacements to the relative displacement of the tip bond. Applying equal and opposite unit loads at its endpoints gives the relative compliance c = h transpose times inverse(K) times h. Cutting the bond changes the matrix by the rank-one term -k h h transpose. Force balance on the remaining network gives a relaxation contribution fcrit squared times c divided by twice (1 - kc). Adding W gives the identity D = W/(1 - kc), a use of the rank-one inverse update [3]. The calculation uses a linear solve for the action of inverse(K). The direct energy subtraction and the compliance calculation require different numerical operations, although they share the assembled network.

## Predictors and finite-size check

The first control predicts W. The second predicts W times A(1), where A = D/W and the argument is patch contrast. Each geometry has its own homogeneous calibration. This makes the control stronger than a single universal multiplier. It knows the notch and specimen sizes as well as the actual stiffness of the breaking spring. Its only missing input is how the surrounding network redistributes force when the patch changes stiffness.

The third predictor uses the local slope of log A against log r at r = 1. Two solves at r = exp(-0.01) and exp(0.01) estimate that slope by a centred difference. The predicted multiplier is A(1) times r raised to the estimated slope. This is a local sensitivity model for a known patch position and geometry. It requires information about which springs are affected. The study tests contrast values that were not used to estimate the slope.

The grid crosses widths 16, 24, 32 and 48 with notch lengths n/4 and n/2, three horizontal offsets, two radii, and seven contrasts from 0.25 to 4. This gives 336 configurations. Forty-eight homogeneous rows serve as controls and 288 rows assess transfer. Some homogeneous rows are physically identical because the nominal patch has unit contrast. The summary gives each design cell equal weight. There is no experimental sampling distribution behind the reported percentages.

## Where does a uniform calibration lose accuracy?

Table 1 gives mean absolute error relative to the directly calculated released energy, in percent. Each contrast averages 48 geometry and patch specifications.

| Contrast | Spring energy alone | Uniform multiplier | Local slope correction |
|---:|---:|---:|---:|
| 0.25 | 55.73 | 6.97 | 1.69 |
| 0.50 | 56.62 | 3.88 | 0.37 |
| 0.80 | 57.34 | 1.31 | 0.03 |
| 1.00 | 57.69 | 0.00 | 0.00 |
| 1.25 | 58.04 | 1.32 | 0.02 |
| 2.00 | 58.73 | 4.05 | 0.15 |
| 4.00 | 59.51 | 7.57 | 0.29 |

Uniform calibration captures most of the surrounding network's contribution. It can still overshoot a specific case appreciably. The largest absolute transfer error, 25.11%, occurs for a 48-node-wide lattice with a 12-bond notch, a tip-centred radius of 1.25 and a fourfold reduction in patch stiffness. The local spring contains 2 units of energy when it reaches its breaking force. Direct calculation releases 3.8385 units. Applying the uniform multiplier overestimates relaxation here. The first-order correction cuts the absolute relative error to 5.21%.

Patch position matters even when its radius is unchanged. Across nonunit contrasts, uniform transfer has 13.88% mean error for a tip-centred radius of 1.25. The same radius yields 0.15% four spacings behind the tip and 0.26% four spacings ahead. Enlarging the concentric radius to 2.75 reduces the centred error to 6.47%. When a small patch covers the cut spring but few alternative load paths, local rupture energy and network relaxation change by different factors. A broad patch changes more surrounding paths as well. Local contrast alone cannot encode that distinction.

The largest corrected error is 5.96%, on a width-16 half-notched lattice with the small soft tip patch. The softening and stiffening results are asymmetric. For the strongest softening the mean corrected error is 1.69%, versus 0.29% for the strongest stiffening. An expansion near unity has limited reach in the soft region even with exact patch membership.

## Boundary effects and physical scope

On a homogeneous half-notched specimen the calculated energy release increases from 1.1516 at width 16 to 1.1715, 1.1809 and 1.1900 at widths 24, 32 and 48. By contrast, the corrected predictor's mean error over nonunit contrasts is 0.46% at width 16 and 0.41% at width 48. A small prediction error is possible even when the absolute release still depends on the boundaries. The lattice spacing is held fixed while the specimen grows. This is a finite-size check and cannot establish continuum mesh convergence.

The study measures a single cut in antiplane shear, with linear scalar springs and a prescribed common breaking force. It omits bending, finite stretch, strand contour length, crack branching and dissipation. The homogeneous multiplier is an internal control for this model. Hartquist's published constitutive scaling is tested on a different quantity, intrinsic fracture energy, and an opening-mode physical comparison would need the corresponding strand force-length curves, density and crack-area normalization [1]. In particular, one cannot compare the percentage errors in Table 1 with a measured toughness error.

## Conclusion and reproducibility

A compact tip patch is a useful first specimen for testing transfer because it maximized the error of homogeneous calibration in this grid. A paired control would put the same patch four cells away and measure the work required for repeated crack advances under fixed grips. Recording the failure sequence would also test whether the prescribed tip cut occurs before off-path rupture.

Running `python run.py` regenerates results.csv, including the imposed grip separation at critical load, the direct release, all predictors and the relative residual of the identity. The largest residual was 6.18e-14. All inputs are generated locally. The study used Python 3.11.9, NumPy 1.26.0 and SciPy 1.11.2.

## Research and writing assistance

Vedant Thaker supplied the research brief and source materials. OpenAI Codex assisted with study design, implementation, execution, literature review and drafting. Clever Humanizer suggested edits to the abstract. The final wording was checked against the analysis.

## References

[1] C. Hartquist, S. Wang, Q. Cui, W. Matusik, B. Deng and X. Zhao. Scaling Law for Intrinsic Fracture Energy of Diverse Stretchable Networks. Physical Review X 15, 011002 (2025). https://doi.org/10.1103/PhysRevX.15.011002

[2] M. Nieves, A. Movchan, I. Jones and G. Mishuris. Propagation of Slepyan's crack in a non-uniform elastic lattice. Journal of the Mechanics and Physics of Solids (2013). https://arxiv.org/abs/1204.5766

[3] J. Sherman and W. J. Morrison. Adjustment of an Inverse Matrix Corresponding to a Change in One Element of a Given Matrix. Annals of Mathematical Statistics 21, 124-127 (1950). https://doi.org/10.1214/aoms/1177729893
