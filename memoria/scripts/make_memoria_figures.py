#!/usr/bin/env python3
"""
Generate PDF-only auxiliary figures and CSV tables for the X17 Typst report.

This script intentionally does not require ROOT.  It reads the supplied
`data/data_pair_creation.txt`, reproduces the simple signal/background sampling
logic used by the Geant4 primary generator, and performs an independent
ray-box geometrical coincidence estimate using the silicon dimensions and the
four geometry layouts defined in `src/DetectorConstruction.cc`.

For the official ROOT tree summaries, run `extract_root_summaries.C` in a local
ROOT-enabled environment.  The figures and tables created here are nevertheless
useful because they are deterministic and document the generator-level model.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

ME = 0.51099895  # MeV/c^2

# Geometry constants in mm, matching DetectorConstruction.cc
SCINT_HALF_THICKNESS = 10.0
SCINT_HALF_WIDTH = 41.0
SCINT_HALF_HEIGHT = 43.0
SI_HALF_THICKNESS = 0.25
SI_HALF_WIDTH = 41.0
SI_HALF_HEIGHT = 43.0
GAP_SI_SCINT = 2.0
SI_HALF_BOX = np.array([SI_HALF_THICKNESS, SI_HALF_WIDTH, SI_HALF_HEIGHT], dtype=float)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def mass_ee_from_kinetic(k1: np.ndarray, k2: np.ndarray, theta_deg: np.ndarray) -> np.ndarray:
    e1 = k1 + ME
    e2 = k2 + ME
    p1 = np.sqrt(np.maximum(0.0, e1 * e1 - ME * ME))
    p2 = np.sqrt(np.maximum(0.0, e2 * e2 - ME * ME))
    c12 = np.cos(np.deg2rad(theta_deg))
    m2 = 2.0 * ME * ME + 2.0 * (e1 * e2 - p1 * p2 * c12)
    return np.sqrt(np.maximum(0.0, m2))


def set_paper_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 0.8,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.figsize": (5.2, 3.35),
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


def hist_overlay(data_series: List[Tuple[np.ndarray, str]], bins: np.ndarray, xlabel: str, ylabel: str,
                 title: str, output: Path, density: bool = True) -> None:
    fig, ax = plt.subplots()
    for values, label in data_series:
        ax.hist(values, bins=bins, histtype="step", linewidth=1.6, density=density, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.22, linewidth=0.5)
    ax.legend(frameon=False)
    fig.savefig(output)
    plt.close(fig)


def direction_from_polar_azimuth(theta_deg: np.ndarray, phi: np.ndarray) -> np.ndarray:
    theta = np.deg2rad(theta_deg)
    return np.column_stack((
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
    ))


def random_isotropic(rng: np.random.Generator, n: int) -> np.ndarray:
    cos_theta = 2.0 * rng.random(n) - 1.0
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - cos_theta * cos_theta))
    phi = 2.0 * np.pi * rng.random(n)
    return np.column_stack((sin_theta * np.cos(phi), sin_theta * np.sin(phi), cos_theta))


def direction_at_opening(axis: np.ndarray, opening_rad: float, rotation_rad: float) -> np.ndarray:
    u = axis / np.linalg.norm(axis)
    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(u @ ref)) > 0.95:
        ref = np.array([1.0, 0.0, 0.0])
    v = np.cross(u, ref)
    v = v / np.linalg.norm(v)
    w = np.cross(v, u)
    w = w / np.linalg.norm(w)
    transverse = v * math.cos(rotation_rad) + w * math.sin(rotation_rad)
    return u * math.cos(opening_rad) + transverse * math.sin(opening_rad)


@dataclass(frozen=True)
class Sample:
    theta_ee: np.ndarray
    theta_e_minus: np.ndarray
    theta_e_plus: np.ndarray
    energy_e_minus: np.ndarray
    energy_e_plus: np.ndarray
    direction_e_minus: np.ndarray
    direction_e_plus: np.ndarray

    @property
    def energy_sum(self) -> np.ndarray:
        return self.energy_e_minus + self.energy_e_plus

    @property
    def mass_ee(self) -> np.ndarray:
        return mass_ee_from_kinetic(self.energy_e_minus, self.energy_e_plus, self.theta_ee)


def generate_signal(table: np.ndarray, n_events: int, seed: int) -> Sample:
    rng = np.random.default_rng(seed)
    rows = table[np.arange(n_events) % len(table)]
    phi0 = 2.0 * np.pi * rng.random(n_events)
    direction_e_minus = direction_from_polar_azimuth(rows[:, 1], phi0)
    direction_e_plus = direction_from_polar_azimuth(rows[:, 3], phi0 + np.pi)
    return Sample(
        theta_ee=rows[:, 0],
        theta_e_minus=rows[:, 1],
        theta_e_plus=rows[:, 3],
        energy_e_minus=rows[:, 2],
        energy_e_plus=rows[:, 4],
        direction_e_minus=direction_e_minus,
        direction_e_plus=direction_e_plus,
    )


def generate_background(n_events: int, seed: int) -> Sample:
    rng = np.random.default_rng(seed)
    theta_scale_deg = 40.0
    theta_max_deg = 180.0
    norm = 1.0 - math.exp(-theta_max_deg / theta_scale_deg)
    theta_ee = -theta_scale_deg * np.log(1.0 - rng.random(n_events) * norm)

    total_kinetic = 17.5
    min_kinetic = 0.20
    available = total_kinetic - 2.0 * min_kinetic
    energy_e_minus = min_kinetic + available * rng.random(n_events)
    energy_e_plus = total_kinetic - energy_e_minus

    direction_e_minus = random_isotropic(rng, n_events)
    cone_phi = 2.0 * np.pi * rng.random(n_events)
    direction_e_plus = np.array([
        direction_at_opening(axis, math.radians(theta), rot)
        for axis, theta, rot in zip(direction_e_minus, theta_ee, cone_phi)
    ])

    theta_e_minus = np.rad2deg(np.arccos(np.clip(direction_e_minus[:, 2], -1.0, 1.0)))
    theta_e_plus = np.rad2deg(np.arccos(np.clip(direction_e_plus[:, 2], -1.0, 1.0)))

    return Sample(
        theta_ee=theta_ee,
        theta_e_minus=theta_e_minus,
        theta_e_plus=theta_e_plus,
        energy_e_minus=energy_e_minus,
        energy_e_plus=energy_e_plus,
        direction_e_minus=direction_e_minus,
        direction_e_plus=direction_e_plus,
    )


def rotz(deg: float) -> np.ndarray:
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def radial_axes(direction: np.ndarray) -> np.ndarray:
    x_axis = direction / np.linalg.norm(direction)
    z_ref = np.array([0.0, 0.0, 1.0])
    y_ref = np.array([0.0, 1.0, 0.0])
    ref = z_ref if abs(float(x_axis @ z_ref)) < 0.95 else y_ref
    y_axis = np.cross(ref, x_axis)
    y_axis = y_axis / np.linalg.norm(y_axis)
    z_axis = np.cross(x_axis, y_axis)
    z_axis = z_axis / np.linalg.norm(z_axis)
    return np.column_stack((x_axis, y_axis, z_axis))


def pad_axes(sign: float = 1.0) -> np.ndarray:
    # local X is the detector normal; local Y/Z span the detector face.
    x_axis = np.array([0.0, 0.0, sign])
    y_axis = np.array([1.0, 0.0, 0.0])
    z_axis = np.array([0.0, sign, 0.0])
    return np.column_stack((x_axis, y_axis, z_axis))


Detector = Tuple[np.ndarray, np.ndarray]


def current_detectors() -> List[Detector]:
    dets: List[Detector] = []
    hex_side = 2.0 * SCINT_HALF_WIDTH
    hex_apothem = 0.5 * math.sqrt(3.0) * hex_side
    scint_radius = 1.03 * (hex_apothem + SCINT_HALF_THICKNESS)
    silicon_radius = scint_radius - SCINT_HALF_THICKNESS - GAP_SI_SCINT - SI_HALF_THICKNESS
    rotations = {0: 0.0, 1: 120.0, 2: 60.0, 3: 0.0, 4: 120.0, 5: 60.0}
    for i in range(6):
        phi = math.radians(i * 360.0 / 6.0)
        center = 1.1 * silicon_radius * np.array([math.cos(phi), math.sin(phi), 0.0])
        dets.append((center, rotz(rotations[i])))
    return dets


def hemisphere_detectors(two_pi: bool) -> List[Detector]:
    dets: List[Detector] = []
    rings = [(28.0, 6, 0.0), (55.0, 10, 18.0), (80.0, 14, 0.0)]
    hemispheres = [True] if two_pi else [True, False]
    for downstream in hemispheres:
        for theta_from_axis_deg, n_phi, phi_offset_deg in rings:
            theta_deg = theta_from_axis_deg if downstream else 180.0 - theta_from_axis_deg
            theta = math.radians(theta_deg)
            for i in range(n_phi):
                phi_deg = phi_offset_deg + i * 360.0 / n_phi
                phi = math.radians(phi_deg)
                direction = np.array([
                    math.sin(theta) * math.cos(phi),
                    math.sin(theta) * math.sin(phi),
                    math.cos(theta),
                ])
                center = 220.0 * direction
                dets.append((center, radial_axes(direction)))
    return dets


def padplane_detectors() -> List[Detector]:
    dets: List[Detector] = []
    nx = ny = 5
    pitch = 90.0
    z_silicon = 180.0
    axes = pad_axes(+1.0)
    for ix in range(nx):
        for iy in range(ny):
            x = (ix - 0.5 * (nx - 1)) * pitch
            y = (iy - 0.5 * (ny - 1)) * pitch
            if abs(x) < 1e-12 and abs(y) < 1e-12:
                continue
            dets.append((np.array([x, y, z_silicon]), axes))
    return dets


def ray_box_hit(direction: np.ndarray, center: np.ndarray, axes: np.ndarray,
                half_box: np.ndarray = SI_HALF_BOX) -> bool:
    # Transform the ray p(t)=t*direction into detector-local coordinates.
    ro = axes.T @ (-center)
    rd = axes.T @ direction
    t_min = -np.inf
    t_max = np.inf
    for j in range(3):
        if abs(rd[j]) < 1e-14:
            if ro[j] < -half_box[j] or ro[j] > half_box[j]:
                return False
            continue
        t1 = (-half_box[j] - ro[j]) / rd[j]
        t2 = (half_box[j] - ro[j]) / rd[j]
        if t1 > t2:
            t1, t2 = t2, t1
        t_min = max(t_min, t1)
        t_max = min(t_max, t2)
        if t_min > t_max:
            return False
    return bool(t_max > 1e-9 and t_min > 0.0)


def hit_mask(directions: np.ndarray, detectors: List[Detector]) -> np.ndarray:
    """Vectorized ray-box intersection for many directions and many silicon boxes."""
    n = directions.shape[0]
    hit = np.zeros(n, dtype=bool)
    for center, axes in detectors:
        # local ray origin relative to the detector and local ray direction.
        ro = axes.T @ (-center)
        rd = directions @ axes

        t_min = np.full(n, -np.inf)
        t_max = np.full(n, np.inf)
        valid = np.ones(n, dtype=bool)

        for j in range(3):
            rj = rd[:, j]
            parallel = np.abs(rj) < 1e-14
            valid &= (~parallel) | ((ro[j] >= -SI_HALF_BOX[j]) & (ro[j] <= SI_HALF_BOX[j]))

            nonparallel = ~parallel
            t1 = np.empty(n)
            t2 = np.empty(n)
            t1[nonparallel] = (-SI_HALF_BOX[j] - ro[j]) / rj[nonparallel]
            t2[nonparallel] = ( SI_HALF_BOX[j] - ro[j]) / rj[nonparallel]
            t1[parallel] = -np.inf
            t2[parallel] = np.inf
            lo = np.minimum(t1, t2)
            hi = np.maximum(t1, t2)
            t_min = np.maximum(t_min, lo)
            t_max = np.minimum(t_max, hi)

        hit |= valid & (t_min <= t_max) & (t_max > 1e-9) & (t_min > 0.0)
    return hit


def coincidence_efficiency(sample: Sample, detectors: List[Detector]) -> Tuple[int, int, int, int, float, float]:
    hit_e_minus = hit_mask(sample.direction_e_minus, detectors)
    hit_e_plus = hit_mask(sample.direction_e_plus, detectors)
    coincidence = hit_e_minus & hit_e_plus
    n = len(coincidence)
    eff = float(np.mean(coincidence)) if n else 0.0
    err = math.sqrt(eff * (1.0 - eff) / n) if n else 0.0
    return int(hit_e_minus.sum()), int(hit_e_plus.sum()), int(coincidence.sum()), n, eff, err


def write_sampling_summary(table: np.ndarray, signal: Sample, background: Sample, output: Path) -> None:
    rows = [
        ("signal", len(signal.theta_ee), signal.theta_ee.mean(), signal.theta_ee.std(ddof=1),
         signal.energy_sum.mean(), signal.energy_sum.std(ddof=1), signal.mass_ee.mean(), signal.mass_ee.std(ddof=1)),
        ("background", len(background.theta_ee), background.theta_ee.mean(), background.theta_ee.std(ddof=1),
         background.energy_sum.mean(), background.energy_sum.std(ddof=1), background.mass_ee.mean(), background.mass_ee.std(ddof=1)),
    ]
    with output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "n", "thetaEE_mean_deg", "thetaEE_std_deg", "energy_sum_mean_MeV",
                         "energy_sum_std_MeV", "massEE_mean_MeV", "massEE_std_MeV"])
        for row in rows:
            writer.writerow(row)


def write_efficiency_tables(signal: Sample, background: Sample, tables_dir: Path, fig_dir: Path) -> None:
    geometries: Dict[str, List[Detector]] = {
        "current": current_detectors(),
        "2pi": hemisphere_detectors(two_pi=True),
        "4pi": hemisphere_detectors(two_pi=False),
        "padplane": padplane_detectors(),
    }

    eff_rows = []
    reco_rows = []
    for geom, detectors in geometries.items():
        for label, sample in (("signal", signal), ("background", background)):
            n_em, n_ep, n_coin, n_gen, eff, err = coincidence_efficiency(sample, detectors)
            eff_rows.append({
                "geometry": geom,
                "sample": label,
                "nGenerated": n_gen,
                "nDetectedElectron": n_em,
                "nDetectedPositron": n_ep,
                "nDetectedCoincidence": n_coin,
                "effCoincidence": eff,
                "effCoincidenceErr": err,
            })

        s_eff = next(r for r in eff_rows if r["geometry"] == geom and r["sample"] == "signal")
        b_eff = next(r for r in eff_rows if r["geometry"] == geom and r["sample"] == "background")
        # Analysis-level detector effects used by apply_detector_effects.C: 90% per particle.
        # Background additionally loses events when one particle falls below the 1 MeV threshold.
        # For the signal table this is evaluated directly from the generated energies.
        sig_threshold = float(np.mean((signal.energy_e_minus >= 1.0) & (signal.energy_e_plus >= 1.0)))
        bkg_threshold = float(np.mean((background.energy_e_minus >= 1.0) & (background.energy_e_plus >= 1.0)))
        sig_reco = s_eff["nDetectedCoincidence"] * 0.90 * 0.90 * sig_threshold
        bkg_reco = b_eff["nDetectedCoincidence"] * 0.90 * 0.90 * bkg_threshold
        total = sig_reco + bkg_reco
        reco_rows.append({
            "geometry": geom,
            "x17_reconstructed_equal_samples": sig_reco,
            "background_reconstructed_equal_samples": bkg_reco,
            "total_reconstructed_equal_samples": total,
            "true_x17_fraction_equal_samples": sig_reco / total if total > 0 else 0.0,
            "background_fraction_equal_samples": bkg_reco / total if total > 0 else 0.0,
        })

    with (tables_dir / "geometrical_efficiencies.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(eff_rows[0].keys()))
        writer.writeheader()
        writer.writerows(eff_rows)

    with (tables_dir / "reconstructed_yield_estimates.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(reco_rows[0].keys()))
        writer.writeheader()
        writer.writerows(reco_rows)

    # Typst fragments keep the main paper readable and avoid manual transcription errors.
    with (tables_dir / "geometrical_efficiencies.typ").open("w") as f:
        f.write("#table(\n")
        f.write("  columns: (1.15fr, 1.20fr, 1.20fr),\n")
        f.write("  inset: 5pt,\n")
        f.write("  align: (left, right, right),\n")
        f.write("  [*Geometry*], [*Signal*], [*Background*],\n")
        for geom in ["current", "2pi", "4pi", "padplane"]:
            s = next(r for r in eff_rows if r["geometry"] == geom and r["sample"] == "signal")
            b = next(r for r in eff_rows if r["geometry"] == geom and r["sample"] == "background")
            f.write(f"  [{geom}], [{100*s['effCoincidence']:.2f} ± {100*s['effCoincidenceErr']:.2f}%], ")
            f.write(f"[{100*b['effCoincidence']:.2f} ± {100*b['effCoincidenceErr']:.2f}%],\n")
        f.write(")\n")

    with (tables_dir / "reconstructed_yields.typ").open("w") as f:
        f.write("#table(\n")
        f.write("  columns: (1.05fr, 1.10fr, 1.10fr, 1.05fr),\n")
        f.write("  inset: 4.5pt,\n")
        f.write("  align: (left, right, right, right),\n")
        f.write("  [*Geometry*], [*$N_X^reco$*], [*$N_b^reco$*], [*true X17 fraction*],\n")
        for row in reco_rows:
            f.write(f"  [{row['geometry']}], [{row['x17_reconstructed_equal_samples']:.0f}], ")
            f.write(f"[{row['background_reconstructed_equal_samples']:.0f}], ")
            f.write(f"[{100*row['true_x17_fraction_equal_samples']:.1f}%],\n")
        f.write(")\n")

    # Bar chart comparing coincidence efficiencies.
    geoms = ["current", "2pi", "4pi", "padplane"]
    x = np.arange(len(geoms))
    width = 0.35
    signal_eff = np.array([next(r for r in eff_rows if r["geometry"] == g and r["sample"] == "signal")["effCoincidence"] for g in geoms])
    bkg_eff = np.array([next(r for r in eff_rows if r["geometry"] == g and r["sample"] == "background")["effCoincidence"] for g in geoms])
    signal_err = np.array([next(r for r in eff_rows if r["geometry"] == g and r["sample"] == "signal")["effCoincidenceErr"] for g in geoms])
    bkg_err = np.array([next(r for r in eff_rows if r["geometry"] == g and r["sample"] == "background")["effCoincidenceErr"] for g in geoms])

    fig, ax = plt.subplots(figsize=(5.4, 3.25))
    ax.bar(x - width / 2, 100.0 * signal_eff, width, yerr=100.0 * signal_err, capsize=2.5, label="X17-like signal")
    ax.bar(x + width / 2, 100.0 * bkg_eff, width, yerr=100.0 * bkg_err, capsize=2.5, label="IPC-like background")
    ax.set_xticks(x, geoms)
    ax.set_ylabel("geometrical coincidence efficiency [%]")
    ax.set_title("Coincidence acceptance by geometry")
    ax.grid(axis="y", alpha=0.22, linewidth=0.5)
    ax.legend(frameon=False)
    fig.savefig(fig_dir / "geometrical_efficiencies.pdf")
    plt.close(fig)


def make_figures(signal: Sample, background: Sample, fig_dir: Path) -> None:
    ensure_dir(fig_dir)
    bins_theta = np.linspace(0, 180, 91)
    bins_energy = np.linspace(0, 20, 81)
    bins_sum = np.linspace(0, 22, 89)
    bins_mass = np.linspace(0, 22, 89)

    hist_overlay([
        (signal.theta_ee, r"$\theta_{ee}$"),
        (signal.theta_e_minus, r"$\theta_{e^-}$"),
        (signal.theta_e_plus, r"$\theta_{e^+}$"),
    ], bins_theta, r"angle [deg]", "normalized events", "Signal sampled angular variables",
        fig_dir / "signal_angles_sampled.pdf")

    hist_overlay([
        (signal.energy_sum, r"$T_{e^-}+T_{e^+}$"),
        (signal.energy_e_minus, r"$T_{e^-}$"),
        (signal.energy_e_plus, r"$T_{e^+}$"),
    ], bins_energy, r"kinetic energy [MeV]", "normalized events", "Signal sampled energy variables",
        fig_dir / "signal_energies_sampled.pdf")

    hist_overlay([
        (background.theta_ee, r"$\theta_{ee}$"),
        (background.theta_e_minus, r"$\theta_{e^-}$"),
        (background.theta_e_plus, r"$\theta_{e^+}$"),
    ], bins_theta, r"angle [deg]", "normalized events", "IPC-like background sampled angular variables",
        fig_dir / "background_angles_sampled.pdf")

    hist_overlay([
        (background.energy_sum, r"$T_{e^-}+T_{e^+}$"),
        (background.energy_e_minus, r"$T_{e^-}$"),
        (background.energy_e_plus, r"$T_{e^+}$"),
    ], bins_energy, r"kinetic energy [MeV]", "normalized events", "IPC-like background sampled energy variables",
        fig_dir / "background_energies_sampled.pdf")

    hist_overlay([
        (signal.theta_ee, "X17-like signal"),
        (background.theta_ee, "IPC-like background"),
    ], bins_theta, r"$\theta_{ee}$ [deg]", "normalized events", "Generated opening-angle templates",
        fig_dir / "signal_background_thetaee_sampled.pdf")

    hist_overlay([
        (signal.mass_ee, "X17-like signal"),
        (background.mass_ee, "IPC-like background"),
    ], bins_mass, r"$m_{ee}$ [MeV/$c^2$]", "normalized events", "Generated invariant-mass templates",
        fig_dir / "signal_background_mass_sampled.pdf")

    hist_overlay([
        (signal.energy_sum, "X17-like signal"),
        (background.energy_sum, "IPC-like background"),
    ], bins_sum, r"$T_{e^-}+T_{e^+}$ [MeV]", "normalized events", "Generated summed kinetic energy",
        fig_dir / "signal_background_energy_sum_sampled.pdf")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--n-events", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()

    project = args.project_root.resolve()
    memoria = project / "memoria"
    figures_sampling = memoria / "figures" / "sampling"
    tables_dir = memoria / "tables"
    ensure_dir(figures_sampling)
    ensure_dir(tables_dir)

    set_paper_style()
    table = np.loadtxt(project / "data" / "data_pair_creation.txt")
    signal = generate_signal(table, args.n_events, args.seed)
    background = generate_background(args.n_events, args.seed)

    make_figures(signal, background, figures_sampling)
    write_sampling_summary(table, signal, background, tables_dir / "sampling_summary.csv")
    write_efficiency_tables(signal, background, tables_dir, figures_sampling)

    print(f"Wrote PDF figures to {figures_sampling}")
    print(f"Wrote CSV/Typst tables to {tables_dir}")


if __name__ == "__main__":
    main()
