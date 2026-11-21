# Preparation record

When one reference solve can replace a new solve for a cut spring network

Vedant Thaker · 21 November 2026

## Original work

Vedant Thaker is the sole author. He designed the original scalar lattice fracture study, wrote its manuscript, developed its original software, and requested the comparisons and checks reported in this revision.

The author reports that the research project began on 2 August 2026. The previous versions and the new revision represent different stages of the work.

## Revision activities

OpenAI Codex assisted with the rank-one and energy-bound calculations, held-out field generation, direct-solve checks, screening and fallback timing comparisons, and reference verification, as well as literature checks and preparation of the revised text. The author set the research question and requested the revision.

## Computational checks

All 2,160 releases lay within their independently computed bounds to numerical precision. In 30 sentinel cases the direct dipole release agreed with a separately assembled before/after energy calculation to within 1.27 × 10⁻¹³ relative error. A verifier checks case counts, bound ordering, the five acceptance counts and timing denominators. The previous independent reconstruction of the circular grid also agrees with the pinned release, grip and predictor fields to within 2.6 × 10⁻¹² relative discrepancy across all 336 rows. These checks address the discrete model.

## Data and version records

The baseline repository contains the circular-patch reference code and 336-row dataset. The new scripts and machine-readable outputs accompany this paper in the revision supplement. Use transfer.py for the calibration and mixed timings, stress_timing.py for the all-fallback timing experiment, and verify_transfer.py for computational checks. The README and pinned requirements specify the invocation and environment. All fields are generated from recorded seeds.

[Original pinned software](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/f943322c44df01ecfcb4f81c1b686315a79a9e6a) · [Revision code and data](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/main/revision/supplement)

The reported checks are computational. No independent personal rerun or line-by-line review of these revisions by the author is asserted here.
