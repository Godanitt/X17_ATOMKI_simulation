import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Optional quick-check script. For final paper figures, use ROOT:
#   bash scripts/make_root_plots.sh

project_dir = Path(__file__).resolve().parents[1]
data_path = project_dir / "data" / "data_pair_creation.txt"
out_dir = project_dir / "plots"
out_dir.mkdir(exist_ok=True)

if not data_path.exists():
    raise FileNotFoundError(f"Input file not found: {data_path}")

data = np.loadtxt(data_path)

theta_ee = data[:, 0]
theta_em = data[:, 1]
theta_ep = data[:, 3]

plt.figure(figsize=(8, 5))
plt.hist(theta_ee, bins=90, range=(0, 180), histtype="step", linewidth=2, label=r"$\theta_{ee}$")
plt.hist(theta_em, bins=90, range=(0, 180), histtype="step", linewidth=1.5, label=r"$\theta_{e^-}$")
plt.hist(theta_ep, bins=90, range=(0, 180), histtype="step", linewidth=1.5, label=r"$\theta_{e^+}$")
plt.xlabel("Angle [deg]")
plt.ylabel("Entries / 2 deg")
plt.xlim(0, 180)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

out = out_dir / "angular_distributions_matplotlib.pdf"
plt.savefig(out)
plt.show()

print(f"Rows read: {len(data)}")
print(f"theta_ee: {theta_ee.min():.2f} deg - {theta_ee.max():.2f} deg")
print(f"PDF saved to {out}")
