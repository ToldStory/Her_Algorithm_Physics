"""Test the small-area curvature scaling limit for SU(2) plaquettes.

Run from repo root:

    python experiments/06_su2_curvature_scaling_limit.py

This experiment is the non-Abelian companion to the U(1) curvature scaling
layer.  It samples a smooth synthetic SU(2) connection, computes ordinary and
Her-reconstructed plaquette holonomies, takes principal SU(2) logarithm vectors, and compares log(W)/h^2 to the analytic non-Abelian curvature vector and its conjugacy-invariant norm.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from her_reanchoring.su2_curvature_scaling import (  # noqa: E402
    estimate_loglog_slope,
    su2_curvature_scaling_table,
)


def main() -> None:
    rows = su2_curvature_scaling_table((6, 8, 10, 12, 16, 20, 24))
    vector_slope = estimate_loglog_slope(
        [row.h for row in rows],
        [row.rms_vector_curvature_error for row in rows],
    )
    norm_slope = estimate_loglog_slope(
        [row.h for row in rows],
        [row.rms_norm_curvature_error for row in rows],
    )

    out_tables = ROOT / "outputs" / "tables"
    out_figures = ROOT / "outputs" / "figures"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figures.mkdir(parents=True, exist_ok=True)

    csv_path = out_tables / "su2_curvature_scaling_limit_summary.csv"
    fieldnames = [
        "n_side",
        "h",
        "n_plaquettes",
        "rms_vector_curvature_error",
        "max_vector_curvature_error",
        "rms_norm_curvature_error",
        "max_norm_curvature_error",
        "max_her_reconstruction_error",
        "max_plaquette_log_norm",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: getattr(row, name) for name in fieldnames})

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(
        [row.h for row in rows],
        [row.rms_vector_curvature_error for row in rows],
        marker="o",
        label="vector error ||log(W)/h^2 - F||",
    )
    ax.loglog(
        [row.h for row in rows],
        [row.rms_norm_curvature_error for row in rows],
        marker="s",
        label="norm error | ||log(W)/h^2|| - ||F|| |",
    )
    ax.set_title("SU(2) curvature scaling from Her-reconstructed plaquettes")
    ax.set_xlabel("lattice spacing h")
    ax.set_ylabel("RMS curvature error")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig_path = out_figures / "su2_curvature_scaling_limit.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    print("SU(2) curvature scaling limit experiment")
    print("case:                         smooth synthetic SU(2) connection")
    print("diagnostic:                   principal log(W)/h^2 vs analytic non-Abelian F_xy")
    print(f"estimated vector log-log slope: {vector_slope:.3f}")
    print(f"estimated norm log-log slope:   {norm_slope:.3f}")
    for row in rows:
        print(
            f"n={row.n_side:2d}, h={row.h:.5f}, "
            f"rms vector err={row.rms_vector_curvature_error:.3e}, "
            f"rms norm err={row.rms_norm_curvature_error:.3e}, "
            f"max Her recon err={row.max_her_reconstruction_error:.3e}, "
            f"max |log W|={row.max_plaquette_log_norm:.3e}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    main()
