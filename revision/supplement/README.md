# Reusing one reference solve for a prescribed cut

## Reproduce

This directory has two experiments. For the November 21 benchmark, use CPython 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. The supplied revision lock has Linux x86_64 wheel hashes for this environment.

```sh
python -m pip install --require-hashes -r revision-cp312-linux-x86_64.lock
python transfer.py
python stress_timing.py
python verify_transfer.py
```

The scripts require no network after installation. Use a separate working copy if you wish to retain the measured times shipped here. The two timing scripts regenerate their wall-clock CSVs, and timings vary across runs. All three numerical-library thread-count variables default to one in these scripts. The reported timings were measured on an AMD EPYC 7763 under the platform recorded in environment.json. Files do not log any personal identifiers.

- held_out.csv has 2,160 equally weighted configurations. The key is (n, notch, family, seed, ratio, sign). Columns include exact released energy, an alternative-path uniform estimate, a first-order log-sensitivity estimate, the center and both ends of the variational bracket, the tip stiffness, and a per-query relative-error certificate (a fraction, not a percentage).
- energy_checks.csv contains 30 paired energy-difference checks (seed 0, r=4, positive sign, one per family and topology).
- timings_ms.csv has five paired repeat rows per topology for an identical randomly selected 90-field mixed sample. Times are milliseconds per query; setup columns record the first calibration and median of nine repeated setups, in milliseconds per calibration. The direct method reuses its structural assembly indices but factors each new stiffness matrix once.
- timings_high_contrast_ms.csv pairs direct and hybrid times for 48 fields per topology, all binary and at amplitude 4. No certificate passes 2% on that panel.

For the archived circular-patch experiment, run `python run.py` and `python verify.py` to regenerate and confirm results.csv (336 configurations). This independent implementation includes competing-bond forces. To compare it with the pinned [public baseline](https://github.com/VEDANT1602/lattice-fracture-calibration/tree/f943322c44df01ecfcb4f81c1b686315a79a9e6a), run `python compare_archive.py /path/to/baseline/results.csv`. The public baseline uses its own requirements and CPython 3.11 environments; the baseline lock identifies those versioned Linux wheels. No DOI or public release is asserted for the newer supplemental experiment.
