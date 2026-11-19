# Estimating the energy released by one spring cut

A cut releases the strain energy in a loaded spring and lets neighboring springs relax. How much relaxation can a calibration on a uniform lattice account for when the crack tip actually sits in a soft patch?

One configuration gives the size of the error. For a 16-wide lattice with an eight-spring notch and a quarter-stiffness patch centered on the tip, the cut spring stores 2.000 energy units. Resolving equilibrium after cutting gives 3.755 units released. The uniform-lattice correction predicts 4.606. Including the calibration slope with respect to patch stiffness gives 3.531. That is the largest relative error from the slope correction in this grid, 5.96%. [Configuration and both predictions](results.csv?plain=1#L58)

## Calibration across the grid

The script loads a square antiplane lattice to unit force in its prescribed tip spring. It holds grip displacement fixed through the cut and checks the change in total energy against an independent rank-one stiffness calculation. The grid covers 336 configurations, including 48 homogeneous controls. [Energy identity](run.py#L58)

| Estimate | Mean absolute relative error over 288 nonuniform configurations |
| --- | ---: |
| Tip spring energy | 57.66% |
| Geometry-matched uniform calibration | 4.18% |
| Calibration with local stiffness slope | 0.43% |

The tabulated quantity is the release from that single prescribed cut. Predicting which spring fails next, or estimating a physical toughness, would require a crack-growth and failure model. [Full results](results.csv) · [Derivation](manuscript.md) · [Paper PDF](manuscript.pdf)

## Reproduce the calculations

Run from the repository root with Python 3.11. The study used Python 3.11.9, NumPy 1.26.0 and SciPy 1.11.2.

```bash
python -m pip install -r requirements.txt
python run.py
python verify.py
```

`run.py` rebuilds the lattices, solves the equilibrium systems and writes `results.csv`. `verify.py` uses the Python standard library to check the configuration counts, recompute the error table and verify agreement between the two energy calculations.

Repository maintained by **Vedant Thaker**.
