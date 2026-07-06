"""First U(1) gauge re-anchoring experiment.

Run from repo root:

    python experiments/02_u1_lattice_gauge.py

Outputs:
    outputs/tables/u1_gauge_reanchoring_summary.csv
    outputs/figures/u1_perturbed_residual_histogram.png
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from her_reanchoring.u1_gauge import (
    angle_errors,
    gauge_residual_covariance_error,
    gauge_transform,
    perturb_links,
    random_pure_gauge,
    residual_participation_scores,
    residual_summary,
    top_scored_edges,
    triangle_residuals,
)


ROOT = Path(__file__).resolve().parents[1]
OUT_TABLES = ROOT / "outputs" / "tables"
OUT_FIGS = ROOT / "outputs" / "figures"
OUT_TABLES.mkdir(parents=True, exist_ok=True)
OUT_FIGS.mkdir(parents=True, exist_ok=True)


def summary_row(name: str, U: np.ndarray) -> dict[str, str | int | float]:
    s = residual_summary(U)
    return {
        "case": name,
        "n_sites": s.n_sites,
        "n_triangles_ordered": s.n_triangles,
        "max_chord_error": s.max_chord_error,
        "mean_chord_error": s.mean_chord_error,
        "max_angle_error_rad": s.max_angle_error,
        "mean_angle_error_rad": s.mean_angle_error,
    }


def main() -> None:
    n = 12
    alpha, U = random_pure_gauge(n, seed=42)

    lambdas = np.random.default_rng(43).uniform(-np.pi, np.pi, size=n)
    U_gauge = gauge_transform(U, lambdas)

    corrupt_edges = [(0, 1), (2, 5), (3, 7), (8, 10)]
    U_bad, deltas = perturb_links(U, corrupt_edges, strength=0.9, seed=44)
    U_bad_gauge = gauge_transform(U_bad, lambdas)

    rows = [
        summary_row("pure_gauge", U),
        summary_row("pure_gauge_after_local_gauge_transform", U_gauge),
        summary_row("perturbed_sparse_links", U_bad),
        summary_row("perturbed_after_local_gauge_transform", U_bad_gauge),
    ]

    cov_err_pure = gauge_residual_covariance_error(U, lambdas)
    cov_err_bad = gauge_residual_covariance_error(U_bad, lambdas)

    csv_path = OUT_TABLES / "u1_gauge_reanchoring_summary.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
        writer.writerow({})
        writer.writerow({"case": "gauge_invariance_error_pure", "max_chord_error": cov_err_pure})
        writer.writerow({"case": "gauge_invariance_error_perturbed", "max_chord_error": cov_err_bad})

    triples, R_bad = triangle_residuals(U_bad)
    err = angle_errors(R_bad)
    fig_path = OUT_FIGS / "u1_perturbed_residual_histogram.png"
    plt.figure(figsize=(7, 4.5))
    plt.hist(err, bins=40)
    plt.xlabel("absolute triangle residual angle |arg(R_ijk)| [radians]")
    plt.ylabel("ordered triangle count")
    plt.title("U(1) Her re-anchoring residuals after sparse link perturbations")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=160)
    plt.close()

    scores = residual_participation_scores(U_bad)
    top = top_scored_edges(scores, len(corrupt_edges) * 2)
    true = {tuple(sorted(e)) for e in corrupt_edges}
    found_top_m = {edge for edge, _ in top[: len(corrupt_edges)]} & true
    found_top_2m = {edge for edge, _ in top} & true

    print("U(1) gauge re-anchoring experiment")
    print("===================================")
    for row in rows:
        print(
            f"{row['case']}: max_angle={row['max_angle_error_rad']:.3e}, "
            f"mean_angle={row['mean_angle_error_rad']:.3e}"
        )
    print(f"gauge invariance error, pure:      {cov_err_pure:.3e}")
    print(f"gauge invariance error, perturbed: {cov_err_bad:.3e}")
    print(f"corrupt edge phase deltas: {deltas}")
    print(f"top-m recall:  {len(found_top_m)}/{len(true)}")
    print(f"top-2m recall: {len(found_top_2m)}/{len(true)}")
    print("top scored edges:")
    for edge, score in top:
        mark = "TRUE" if edge in true else "----"
        print(f"  {edge}: {score:.6f} {mark}")
    print(f"wrote {csv_path.relative_to(ROOT)}")
    print(f"wrote {fig_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
