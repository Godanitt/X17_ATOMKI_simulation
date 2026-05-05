# X17 sampled simulation

This version uses `data/data_pair_creation.txt` as the only kinematic input.
The input columns are interpreted as:

```text
thetaEE_deg  thetaEe_deg  energyEe_MeV  thetaEp_deg  energyEp_MeV
```

The output simulation file is `sampled.root` with exactly two trees:

```text
generated   one row per generated e-/e+ event
hits        one row per particle/volume detector crossing
```

`hits.pdg = 11` is the electron e- and `hits.pdg = -11` is the positron e+.
`hits.volumeID = 0` is the silicon layer and `hits.volumeID = 1` is the scintillator.

Run simulation:

```bash
bash scripts/run_generated.sh 10000
```

Visualize the new detector geometry:

```bash
bash scripts/build.sh
```

Analyze hit coincidences:

```bash
bash scripts/analyze_hits.sh sampled.root analysis_hits.root 0
```

The last argument chooses the hit volume used to reconstruct coincidences:

```text
0 = SiliconStripLV
1 = ScintillatorLV
-1 = first hit in any detector volume
```

The analysis writes `analysis_hits.root`, including a `detected` tree with one row per event where both e- and e+ were detected.
