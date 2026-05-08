# X17 ATOMKI simulation: signal, background and detector-effects chain

This package keeps the Geant4 output deliberately minimal:

```text
sampled.root
├── generated   # one row per generated e-/e+ event
└── hits        # one row per particle/volume crossing
```

The `hits` tree stores both electrons and positrons in the same tree. Use
`pdg == 11` for the electron and `pdg == -11` for the positron.

## Signal generation

The signal generator reads the supplied table directly:

```text
data/data_pair_creation.txt
col1 = thetaEE  [deg]
col2 = thetaEe  [deg]
col3 = energyEe [MeV]
col4 = thetaEp  [deg]
col5 = energyEp [MeV]
```

Run:

```bash
bash scripts/run_generated.sh 10000 sampled.root
bash scripts/analyze_hits.sh sampled.root analysis_hits.root 0
bash scripts/apply_detector_effects.sh analysis_hits.root analysis_hits_detector_effects.root
```

The analysis output is:

```text
analysis_hits.root
└── detected

analysis_hits_detector_effects.root
└── detected_detector_effects
```

## Background generation

The background is a simple smooth IPC-like toy model. It is not tuned to data. It is meant
as a broad, falling opening-angle template that goes through exactly the same geometry,
hit reconstruction and detector-effects steps as the signal.

Run:

```bash
bash scripts/run_background.sh 10000 background_sampled.root
bash scripts/analyze_background_hits.sh background_sampled.root background_analysis_hits.root 0
bash scripts/apply_background_detector_effects.sh background_analysis_hits.root background_analysis_hits_detector_effects.root
```

Internally this uses:

```bash
x17sim --mode background
```

The background sample still writes exactly two Geant4 trees:

```text
background_sampled.root
├── generated
└── hits
```

## Detector effects

Detector effects are applied at analysis level, not by modifying the Geant4 truth file.
By default:

```text
efficiency per particle       = 0.90
sigma theta                   = 2.0 deg
sigma phi                     = 2.0 deg
relative energy resolution    = 0.05
energy threshold              = 1.0 MeV
```

Override them as:

```bash
bash scripts/apply_detector_effects.sh input.root output.root 0.90 2.0 2.0 0.05 1.0 12345
```

The detector-effects tree includes useful reconstructed observables:

```text
thetaEE_reco_deg
energyEE_reco_MeV
massEE_reco_MeV
```

where the invariant mass uses the stored kinetic energies and the electron mass.

## Final signal/background comparison

After producing both signal and background detector-effects files:

```bash
bash scripts/final_signal_background_analysis.sh \
  analysis_hits_detector_effects.root \
  background_analysis_hits_detector_effects.root \
  final_signal_background.root \
  1.0 1.0
```

The last two numbers are arbitrary scale factors for the expected signal and background
yields. The output file contains template histograms for:

```text
thetaEE_reco_deg
energyEE_reco_MeV
massEE_reco_MeV
signal + background
normalized shapes
bin-by-bin S/sqrt(B)
```

## Visualization

To compile and open the geometry viewer:

```bash
bash scripts/build.sh
```
