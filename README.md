# X17 ATOMKI simulation: signal, background, detector effects and geometry scans

This package keeps the Geant4 output deliberately minimal:

```text
sampled.root
├── generated   # one row per generated e-/e+ event
└── hits        # one row per particle/volume crossing
```

The `hits` tree stores both electrons and positrons in the same tree. Use
`pdg == 11` for the electron and `pdg == -11` for the positron.

## Predefined geometries

The detector geometry is selected at runtime. No radius or detector-count parameters have
to be supplied manually. The `2pi` option is a semisphere, not a cylindrical barrel.

Available names:

```text
current   # original ATOMKI-like 6-telescope geometry
2pi       # semispherical downstream 2π coverage (+Z hemisphere), with beam-axis hole
4pi       # two mirrored semispherical shells, approximately 4π coverage, with beam-axis holes
padplane  # downstream pad plane in front of the beam, with central pad removed
```


### Orientation and beam-axis hole

All telescope modules are oriented so that the thin local axis of the silicon and scintillator boxes points radially away from the target. In practice, each detector face is tangent to the corresponding shell or plane, and the scintillator is placed behind the silicon along the same radial normal.

The `2pi`, `4pi` and `padplane` layouts deliberately leave the beam line empty: there is no telescope exactly on the `+Z` or `-Z` axis, and the pad plane has no central pad at `(x, y) = (0, 0)`. This avoids an unphysical detector directly in the incoming beam direction.

Examples. These commands build the executable and open the Geant4 viewer for the selected geometry:

```bash
./build.sh current
./build.sh 2pi
./build.sh 4pi
./build.sh padplane
```

To build without opening the viewer:

```bash
./build.sh 4pi --build-only
```

If CMake complains about a stale build directory after moving/unzipping the project, use:

```bash
./build.sh 4pi --clean
```

The executable also accepts the geometry directly:

```bash
./build/x17sim -m macros/run_sampled.mac \
  -i data/data_pair_creation.txt \
  -o results/2pi/sampled.root \
  --mode signal \
  --geometry 2pi
```

A macro can also select the geometry before `/run/initialize`:

```text
/x17/geometry/select 4pi
/run/initialize
/run/beamOn 10000
```

## Full pipeline over all geometries

Run the full signal/background/detector-effects/fit chain for the four predefined geometries. This is batch mode and does not open the viewer; use `./build.sh <geometry>` to preview a geometry interactively.

```bash
bash scripts/run_all.sh all 10000 0
```

In `all` mode the script explicitly runs:

```text
current  2pi  4pi  padplane
```

It attempts all four geometries and prints a final summary. If one geometry fails,
the remaining geometries are still attempted and the script exits non-zero at the end.

The default also runs all geometries:

```bash
bash scripts/run_all.sh
```

To run only one geometry:

```bash
bash scripts/run_all.sh current 10000 0
bash scripts/run_all.sh 2pi 10000 0
bash scripts/run_all.sh 4pi 10000 0
bash scripts/run_all.sh padplane 10000 0
```

Arguments:

```text
$1 = geometry name or all
$2 = number of events per signal/background sample
$3 = hit volumeID used by analysis/analyze_hits.C
```

`volumeID` convention:

```text
 0 = SiliconStripLV
 1 = ScintillatorLV
-1 = any detector volume
```

Each geometry writes into its own folder:

```text
results/
├── current/
├── 2pi/
├── 4pi/
└── padplane/
```

For each geometry the main files are:

```text
results/<geometry>/sampled.root
results/<geometry>/background_sampled.root
results/<geometry>/analysis_hits.root
results/<geometry>/background_analysis_hits.root
results/<geometry>/analysis_hits_detector_effects.root
results/<geometry>/background_analysis_hits_detector_effects.root
results/<geometry>/final_signal_background.root
results/<geometry>/signal_background_fit.root
results/<geometry>/ipc_excess_fit.root
results/<geometry>/config.txt
```

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

Run a single signal sample with a selected geometry:

```bash
bash scripts/run_generated.sh 10000 results/2pi/sampled.root data/data_pair_creation.txt 2pi results/2pi/logs
```

## Background generation

The background is a simple smooth IPC-like toy model. It is not tuned to data. It is meant
as a broad, falling opening-angle template that goes through exactly the same geometry,
hit reconstruction and detector-effects steps as the signal.

Run a single background sample with a selected geometry:

```bash
bash scripts/run_background.sh 10000 results/2pi/background_sampled.root 2pi results/2pi/logs
```

Internally this uses:

```bash
x17sim --mode background --geometry 2pi
```

## Geometrical efficiency in ROOT files

`analysis/analyze_hits.C` now writes a `summary` tree to each analysis ROOT file:

```text
analysis_hits.root
├── detected
└── summary
```

The `summary` tree contains:

```text
geometry
selectedVolumeID
nGenerated
nHitRows
nDetectedElectron
nDetectedPositron
nDetectedCoincidence
effElectron
effPositron
effCoincidence
```

The key geometrical efficiency for pair reconstruction is:

```text
effCoincidence = nDetectedCoincidence / nGenerated
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

Override them with environment variables before calling `run_all.sh`:

```bash
EFFICIENCY=0.85 SIGMA_THETA=3.0 bash scripts/run_all.sh 4pi 10000 0
```

The detector-effects tree includes useful reconstructed observables:

```text
thetaEE_reco_deg
energyEE_reco_MeV
massEE_reco_MeV
```

where the invariant mass uses the stored kinetic energies and the electron mass.

## Final signal/background comparison

After producing both signal and background detector-effects files, the pipeline creates:

```text
results/<geometry>/final_signal_background.root
results/<geometry>/signal_background_fit.root
results/<geometry>/ipc_excess_fit.root
```

The IPC-like excess fit builds pseudo-data from signal + background and extracts an excess
histogram by subtracting the fitted smooth background. The excess histogram is therefore an
extracted X17-like signal shape, while the pure signal MC remains in:

```text
results/<geometry>/analysis_hits_detector_effects.root
```
