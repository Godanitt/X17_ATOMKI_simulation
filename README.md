# X17 ATOMKI-like Geant4 simulation

This version is arranged around a clean reconstruction workflow:

1. generate an X17-like e-/e+ pair with a known opening angle;
2. validate an ideal reconstruction from the generated directions;
3. add detector effects separately using silicon hit positions and scintillator energy deposits;
4. write ROOT ntuples and make paper-style PDF plots.

## Run everything

From the project root:

```bash
bash scripts/run_all.sh
```

The default is intentionally safe for a laptop:

```txt
EVENTS=10000
THREADS=2
BUILD JOBS=4
```

Increase these only after checking that your machine does not overheat.

## Output ROOT files

```txt
x17_output.root      Geant4 ntuples: events, hits
x17_analysis.root    Offline reconstruction ntuple: reco
```

## Main plots

```txt
plots/thetaee_gen_all.pdf
plots/thetaee_generated_ideal_geometry_detector.pdf
plots/thetaee_gen_reco.pdf
plots/thetaee_ideal_resolution.pdf
plots/thetaee_detector_resolution.pdf
plots/thetaee_resolution_comparison.pdf
plots/thetaee_reco_vs_gen.pdf
plots/acceptance_ideal_geometry_detector_vs_thetaee.pdf
plots/acceptance_vs_thetaee.pdf
plots/scint_energy_em_ep.pdf
plots/scint_energy_sum.pdf
plots/scint_energy_asymmetry.pdf
plots/thetaee_vs_Esum_reco.pdf
plots/thetaee_vs_Y_reco.pdf
plots/cutflow.pdf
plots/hit_detector_occupancy.pdf
plots/hit_volume_breakdown.pdf
```

## What the reconstruction levels mean

```txt
Generated input:
    theta_ee from the input table.

Ideal reconstructed:
    theta_ee reconstructed from the generated unit vectors.
    This should lie exactly on the diagonal and have near-zero residual.

Geometrical accepted:
    events where both leptons hit the silicon planes.

Detector effects:
    geometrically accepted events with scintillator energy for both leptons,
    including 3 mm silicon strip quantization and 5% scintillator energy smearing.
```

