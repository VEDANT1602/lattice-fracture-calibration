# When a stiffness query can stop before a point solve

Numerical benchmark by Vedant Thaker, Independent Researcher, Jeffersonville, Indiana, USA. Revised November 24, 2026.

Whether to factor a cut spring network depends on what the caller needs. On a 42,848-unknown network, classical trial-energy bounds settled 854 of 1,024 decisions about whether the response exceeded the uniform-network value. Those bounds alone reached 2% relative accuracy on none of 48 strong fields in a separate value-estimation panel. For eight additional strong fields, a certified scalar estimate took 244.4 ms in total, including setup, against 283.5 ms for full solves with reused-pattern Cholesky. For the first field alone, Cholesky was faster (71.9 versus 91.9 ms). The scalar method first recovered its setup cost on the fourth field in seed order. The benchmark quantifies the cost of these two requested outputs under a fixed load and prescribed cut.

[Paper (PDF)](manuscript.pdf) · [Editable text](manuscript.md) · [Methods and study files](study/README.md) · [Manuscript history](study/VERSION.md)

[Author-owned software license](LICENSE) · [Third-party terms](LICENSING.md) · [Previous version](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/45bda5e498f6dff7bac53908eb15554cc19d9c74)
