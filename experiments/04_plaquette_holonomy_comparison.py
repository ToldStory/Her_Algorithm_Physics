"""Compare Her triangle residuals with standard plaquette holonomy.

Run from repo root:

    python experiments/04_plaquette_holonomy_comparison.py

The experiment is deliberately conservative.  It verifies exact algebraic
compatibility between square plaquette holonomy and Her residuals on the two
triangles formed by a diagonal split.  It does not claim that Her residuals are
identical to curvature.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from her_reanchoring.plaquette_holonomy import (  # noqa: E402
    her_plaquette_reconstruction_error,
    her_triangle_holonomy_conjugacy_error,
    plaquette_angles,
    plaquette_holonomies,
    plaquette_participation_scores,
    plaquette_summary,
    square_plaquettes,
    su2_plaquette_covariance_error,
    su2_plaquette_trace_invariance_error,
    top_scored_plaquettes,
    triangle_factorization_error,
    u1_plaquette_gauge_invariance_error,
)
from her_reanchoring.u1_gauge import (  # noqa: E402
    perturb_links as perturb_u1_links,
    random_pure_gauge as random_u1_pure_gauge,
)
from her_reanchoring.su2_gauge import (  # noqa: E402
    perturb_links as perturb_su2_links,
    random_pure_gauge as random_su2_pure_gauge,
    random_su2_stack,
)


def max_over_plaquettes(fn, U, plaquettes, group):
    return max(float(fn(U, p, group)) for p in plaquettes)


def main() -> None:
    nx, ny = 4, 4
    n_sites = nx * ny
    plaquettes = square_plaquettes(nx, ny)

    # Sparse perturbations on nearest-neighbor square edges.
    # These affect only a small number of plaquettes, making the histogram
    # interpretable as a defect-localization signal rather than a thermal field.
    perturbed_edges = [(1, 2), (5, 6), (9, 13)]

    _, U1_pure = random_u1_pure_gauge(n_sites, seed=100)
    U1_pert, u1_deltas = perturb_u1_links(
        U1_pure,
        edges=perturbed_edges,
        strength=0.9,
        seed=101,
        reciprocal=True,
    )

    _, U2_pure = random_su2_pure_gauge(n_sites, seed=200)
    U2_pert, su2_metadata = perturb_su2_links(
        U2_pure,
        edges=perturbed_edges,
        strength=0.9,
        seed=201,
        reciprocal=True,
    )

    rng = np.random.default_rng(300)
    lambdas = rng.uniform(-np.pi, np.pi, size=n_sites)
    H = random_su2_stack(n_sites, seed=301)

    u1_pure_summary = plaquette_summary(U1_pure, nx, ny, "u1")
    u1_pert_summary = plaquette_summary(U1_pert, nx, ny, "u1")
    su2_pure_summary = plaquette_summary(U2_pure, nx, ny, "su2")
    su2_pert_summary = plaquette_summary(U2_pert, nx, ny, "su2")

    u1_pure_split_error = max_over_plaquettes(triangle_factorization_error, U1_pure, plaquettes, "u1")
    u1_pert_split_error = max_over_plaquettes(triangle_factorization_error, U1_pert, plaquettes, "u1")
    su2_pure_split_error = max_over_plaquettes(triangle_factorization_error, U2_pure, plaquettes, "su2")
    su2_pert_split_error = max_over_plaquettes(triangle_factorization_error, U2_pert, plaquettes, "su2")

    u1_pure_her_error = max_over_plaquettes(her_plaquette_reconstruction_error, U1_pure, plaquettes, "u1")
    u1_pert_her_error = max_over_plaquettes(her_plaquette_reconstruction_error, U1_pert, plaquettes, "u1")
    su2_pure_her_error = max_over_plaquettes(her_plaquette_reconstruction_error, U2_pure, plaquettes, "su2")
    su2_pert_her_error = max_over_plaquettes(her_plaquette_reconstruction_error, U2_pert, plaquettes, "su2")

    u1_conjugacy_error = max(
        her_triangle_holonomy_conjugacy_error(U1_pert, a, b, c, "u1")
        for a, b, c, _ in plaquettes
    )
    su2_conjugacy_error = max(
        her_triangle_holonomy_conjugacy_error(U2_pert, a, b, c, "su2")
        for a, b, c, _ in plaquettes
    )

    u1_gauge_error = u1_plaquette_gauge_invariance_error(U1_pert, lambdas, nx, ny)
    su2_covariance_error = su2_plaquette_covariance_error(U2_pert, H, nx, ny)
    su2_trace_error = su2_plaquette_trace_invariance_error(U2_pert, H, nx, ny)

    u1_scores = plaquette_participation_scores(U1_pert, nx, ny, "u1")
    su2_scores = plaquette_participation_scores(U2_pert, nx, ny, "su2")

    rows = [
        {
            "case": "u1_pure",
            "n_plaquettes": u1_pure_summary.n_plaquettes,
            "max_plaquette_angle": u1_pure_summary.max_plaquette_angle,
            "mean_plaquette_angle": u1_pure_summary.mean_plaquette_angle,
            "max_triangle_factorization_error": u1_pure_split_error,
            "max_her_reconstruction_error": u1_pure_her_error,
            "gauge_or_covariance_error": "",
            "trace_invariance_error": "",
            "conjugacy_error": "",
        },
        {
            "case": "u1_perturbed",
            "n_plaquettes": u1_pert_summary.n_plaquettes,
            "max_plaquette_angle": u1_pert_summary.max_plaquette_angle,
            "mean_plaquette_angle": u1_pert_summary.mean_plaquette_angle,
            "max_triangle_factorization_error": u1_pert_split_error,
            "max_her_reconstruction_error": u1_pert_her_error,
            "gauge_or_covariance_error": u1_gauge_error,
            "trace_invariance_error": "",
            "conjugacy_error": u1_conjugacy_error,
        },
        {
            "case": "su2_pure",
            "n_plaquettes": su2_pure_summary.n_plaquettes,
            "max_plaquette_angle": su2_pure_summary.max_plaquette_angle,
            "mean_plaquette_angle": su2_pure_summary.mean_plaquette_angle,
            "max_triangle_factorization_error": su2_pure_split_error,
            "max_her_reconstruction_error": su2_pure_her_error,
            "gauge_or_covariance_error": "",
            "trace_invariance_error": "",
            "conjugacy_error": "",
        },
        {
            "case": "su2_perturbed",
            "n_plaquettes": su2_pert_summary.n_plaquettes,
            "max_plaquette_angle": su2_pert_summary.max_plaquette_angle,
            "mean_plaquette_angle": su2_pert_summary.mean_plaquette_angle,
            "max_triangle_factorization_error": su2_pert_split_error,
            "max_her_reconstruction_error": su2_pert_her_error,
            "gauge_or_covariance_error": su2_covariance_error,
            "trace_invariance_error": su2_trace_error,
            "conjugacy_error": su2_conjugacy_error,
        },
    ]

    out_tables = ROOT / "outputs" / "tables"
    out_figures = ROOT / "outputs" / "figures"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figures.mkdir(parents=True, exist_ok=True)

    csv_path = out_tables / "plaquette_holonomy_comparison_summary.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    W1 = plaquette_holonomies(U1_pert, plaquettes, "u1")
    W2 = plaquette_holonomies(U2_pert, plaquettes, "su2")
    u1_angles = plaquette_angles(W1, "u1")
    su2_angles = plaquette_angles(W2, "su2")

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(plaquettes))
    width = 0.35
    ax.bar(x - width / 2, u1_angles, width, label="U(1)")
    ax.bar(x + width / 2, su2_angles, width, label="SU(2)")
    ax.set_title("Plaquette holonomy angles after sparse link perturbations")
    ax.set_xlabel("plaquette index")
    ax.set_ylabel("gauge-invariant plaquette angle [rad]")
    ax.set_xticks(x)
    ax.legend()
    fig.tight_layout()
    fig_path = out_figures / "u1_su2_plaquette_angles.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    print("Plaquette / holonomy comparison experiment")
    print(f"grid:                         {nx} x {ny}")
    print(f"n plaquettes:                 {len(plaquettes)}")
    print(f"perturbed edges:              {perturbed_edges}")
    print(f"U(1) perturbations:           {u1_deltas}")
    print(f"SU(2) perturbation metadata:  {su2_metadata}")
    print(f"U(1) pure max angle:          {u1_pure_summary.max_plaquette_angle:.3e}")
    print(f"U(1) perturbed max angle:     {u1_pert_summary.max_plaquette_angle:.3e}")
    print(f"U(1) Her reconstruction err:  {u1_pert_her_error:.3e}")
    print(f"U(1) gauge invariance err:    {u1_gauge_error:.3e}")
    print(f"SU(2) pure max angle:         {su2_pure_summary.max_plaquette_angle:.3e}")
    print(f"SU(2) perturbed max angle:    {su2_pert_summary.max_plaquette_angle:.3e}")
    print(f"SU(2) Her reconstruction err: {su2_pert_her_error:.3e}")
    print(f"SU(2) covariance err:         {su2_covariance_error:.3e}")
    print(f"SU(2) trace invariance err:   {su2_trace_error:.3e}")
    print(f"top U(1) plaquettes:          {top_scored_plaquettes(u1_scores, 4)}")
    print(f"top SU(2) plaquettes:         {top_scored_plaquettes(su2_scores, 4)}")
    print(f"wrote {csv_path}")
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    main()
