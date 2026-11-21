# Choosing a solver for repeated stiffness queries on a cut spring network

Technical report by Vedant Thaker, Independent Researcher, Jeffersonville, Indiana. Correspondence: vedant.thaker1602@gmail.com.

A stiffness-field study may require the value of a local response or just the fraction of fields exceeding a reference. Those requests need different stopping rules. A scalar square network with 42,848 unknowns provides a repeated-query benchmark. When 16 specified springs vary, a reusable Woodbury update gives the direct-solve result in 0.14 ms per field after setup. A trial-field bound resolves a 2% relative-accuracy request for all 48 mildly heterogeneous distributed fields tested and none of 48 strongly heterogeneous fields. Even at strong contrast, the same bound classifies 854 of 1,024 fresh fields relative to a uniform reference without iteration. The location of the changing springs and the requested output both affect solver choice.

- [Read the paper](manuscript.pdf)
- [Read the plain-text manuscript](manuscript.md)
- [Reproduce the analysis and inspect the source and data records](study/README.md)
- [Software license](LICENSE) and [third-party terms](LICENSING.md)
- [Preparation record](study/PREPARATION.md)
- [Previous manuscript and supplements](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/85bf05d1b58dde306f5fca88ef20fbf640c7a7df)
