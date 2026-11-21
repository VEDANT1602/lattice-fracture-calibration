# When a stiffness query can stop before a point solve

Numerical methods report by Vedant Thaker, Independent Researcher, Jeffersonville, Indiana, USA.

**Research period (author-reported)** · August 2–13, 2026 &nbsp; | &nbsp; **Manuscript revised** · November 21, 2026

Estimating a response to a specified accuracy and deciding whether it crosses a threshold can require different work. On a cut spring network with 42,848 free nodes, trial energy bounds settled 854 of 1,024 threshold queries for strongly heterogeneous fields. The same bounds met a 2% point-value tolerance on none of 48 strong fields in a separate panel. Sparse Cholesky with a reusable symbolic factor was a substantially faster direct comparator than fresh SuperLU. On eight other strong fields of this size, median query times were 30.7 ms for reused-pattern Cholesky, 138.8 ms for SuperLU and 20.7 ms for a 2% output certificate. Rebuilding a classical algebraic multigrid hierarchy and converging preconditioned conjugate gradients took 161.5 ms. Five specimen sizes establish the range of this comparison. The practical choice depends on the requested output, the contrast and the number of queries that can share setup.

[Paper (PDF)](manuscript.pdf) · [Editable text](manuscript.md) · [Methods and study files](study/README.md) · [Manuscript history](study/VERSION.md)

[Author-owned software license](LICENSE) · [Third-party terms](LICENSING.md) · [Previous version](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/c4c03fe97b7df961f28efec87e13142cbda0dfb9)
