"""Test the small-area curvature scaling limit for Her-reconstructed plaquettes.

Run from repo root:

    python experiments/05_curvature_scaling_limit.py

This experiment does not claim that Her residuals are curvature.  It checks a
more precise statement: in a smooth synthetic U(1) connection, ordinary
plaquette holonomy has the expected curvature-times-area scaling, and the
Her-reconstructed plaquette holonomy follows the same scaling because it
reconstructs the plaquette holonomy exactly.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from her_reanchoring.curvature_scaling import (  # noqa: E402
    curvature_scaling_table,
    estimate_loglog_slope,
)


def main() -> None:
    rows = curvature_scaling_table((8, 12, 16, 24, 32, 48))
    slope = estimate_loglog_slope(
        [row.h for row in rows],
        [row.rms_curvature_error for row in rows],
    )

    out_tables = ROOT / "outputs" / "tables"
    out_figures = ROOT / "outputs" / "figures"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figures.mkdir(parents=True, exist_ok=True)

    csv_path = out_tables / "curvature_scaling_limit_summary.csv"
    fieldnames = [
        "n_side",
        "h",
        "n_plaquettes",
        "rms_curvature_error",
        "max_abs_curvature_error",
        "rms_angle_error",
        "max_abs_angle_error",
        "max_her_reconstruction_error",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: getattr(row, name) for name in fieldnames})

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog([row.h for row in rows], [row.rms_curvature_error for row in rows], marker="o")
    ax.set_title("U(1) curvature scaling from Her-reconstructed plaquettes")
    ax.set_xlabel("lattice spacing h")
    ax.set_ylabel("RMS error in angle/h^2 versus analytic F_xy")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig_path = out_figures / "curvature_scaling_limit.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    print("Curvature scaling limit experiment")
    print("case:                         smooth synthetic U(1) connection")
    print("diagnostic:                   plaquette angle / h^2 vs analytic F_xy")
    print(f"estimated log-log slope:      {slope:.3f}")
    for row in rows:
        print(
            f"n={row.n_side:2d}, h={row.h:.5f}, "
            f"rms curvature err={row.rms_curvature_error:.3e}, "
            f"max Her recon err={row.max_her_reconstruction_error:.3e}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    main()
