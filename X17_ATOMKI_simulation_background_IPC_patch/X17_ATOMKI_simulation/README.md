# X17 ATOMKI-like Geant4 simulation

This version is arranged around a clean reconstruction workflow:

1. generate a mixed physics sample: dominant IPC-like smooth e-/e+ pairs plus a smaller X17-like signal;
2. validate an ideal reconstruction from the generated directions;
3. add detector effects separately using silicon hit positions, scintillator energy deposits, and a small electronics-noise term;
4. write ROOT ntuples and make paper-style PDF plots, including separated IPC/signal components.

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
SIGNAL_FRACTION=0.10
IPC_THETA_SLOPE_PER_DEG=0.018
```

This means that, by default, 90% of generated pairs are the IPC-like smooth background and 10% are the X17-like signal. Increase event counts only after checking that your machine does not overheat.

You can change the physics mixture without recompiling:

```bash
SIGNAL_FRACTION=0.05 IPC_THETA_SLOPE_PER_DEG=0.020 bash scripts/run_all.sh
```

The Geant4 executable also reads these optional environment variables directly:

```txt
X17_SIGNAL_FRACTION
X17_IPC_THETA_MIN_DEG
X17_IPC_THETA_MAX_DEG
X17_IPC_THETA_SLOPE_PER_DEG
X17_IPC_ENERGY_ASYMMETRY_MAX
X17_IPC_ENERGY_SCALE_SIGMA
```

## Output ROOT files

```txt
x17_output.root      Geant4 ntuples: events, hits
x17_analysis.root    Offline reconstruction ntuple: reco
```

## Main plots

```txt
plots/thetaee_gen_all.pdf
plots/thetaee_components_generated.pdf
plots/thetaee_components_detector.pdf
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
    total theta_ee distribution from the mixed generator.
    componentID = 0: IPC-like smooth background.
    componentID = 1: X17-like signal from data/data_pair_creation.txt.
    componentID = 2 is reserved for accidental coincidences in a later extension.

Ideal reconstructed:
    theta_ee reconstructed from the generated unit vectors.
    This should lie exactly on the diagonal and have near-zero residual.

Geometrical accepted:
    events where both leptons hit the silicon planes.

Detector effects:
    geometrically accepted events with scintillator energy for both leptons,
    including 3 mm silicon strip quantization, 5% scintillator energy smearing,
    and a small 20 keV electronics-noise term. This noise is treated as a
    detector response effect, not as the main physics background.
```



## Background model implemented here

```txt
Primary background:
    IPC-like smooth e-/e+ pairs. The generator samples theta_ee from a
    monotonically decreasing continuum distribution with no resonant bump.

Secondary background:
    accidental coincidences are intentionally not generated yet. The component
    code 2 is reserved so they can be added later without changing the analysis
    interface.

Electronic noise:
    implemented only as a small detector-effect contribution in the offline
    energy smearing. It is not mixed as a separate event class.
```
