# X17 sampled simulation

This version uses `data/data_pair_creation.txt` as the only kinematic input.
The input columns are interpreted as:

```text
thetaEE_deg  thetaEe_deg  energyEe_MeV  thetaEp_deg  energyEp_MeV
```

The simulation output file is `sampled.root` with exactly two trees:

```text
generated   one row per generated e-/e+ event
hits        one row per particle/volume detector crossing
```

`hits.pdg = 11` is the electron e- and `hits.pdg = -11` is the positron e+.
`hits.volumeID = 0` is the silicon layer and `hits.volumeID = 1` is the scintillator.

## Run simulation

```bash
bash scripts/run_generated.sh 10000
```

## Visualize detector geometry

```bash
bash scripts/build.sh
```

## Analyze ideal hit coincidences

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

## Apply analysis-level detector effects

Detector effects are applied after the ideal coincidence analysis, without modifying `sampled.root`.

```bash
bash scripts/apply_detector_effects.sh analysis_hits.root analysis_hits_detector_effects.root
```

Default detector-effect parameters are:

```text
particle efficiency per particle = 0.90
sigma theta                      = 2.0 deg
sigma phi                        = 2.0 deg
relative energy resolution       = 0.05
energy threshold                 = 1.0 MeV
random seed                      = 12345
```

You can override them as:

```bash
bash scripts/apply_detector_effects.sh \
  analysis_hits.root \
  analysis_hits_detector_effects.root \
  0.90 2.0 2.0 0.05 1.0 12345
```

The output file contains:

```text
analysis_hits_detector_effects.root
└── detected_detector_effects
```

with one row per ideal e-/e+ coincidence that survives efficiency, smearing and energy-threshold cuts.
