"""SU(2) non-Abelian gauge re-anchoring experiment.

Run from repo root:

    python experiments/03_su2_nonabelian_covariance.py

Outputs:
    outputs/tables/su2_nonabelian_covariance_summary.csv
    outputs/figures/su2_perturbed_residual_histogram.png
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from her_reanchoring.su2_gauge import (
    gauge_transform,
    perturb_links,
    random_pure_gauge,
    random_su2_stack,
    residual_covariance_error,
    residual_literal_change,
    residual_participation_scores,
    residual_summary,
    residual_trace_invariance_error,
    su2_rotation_angles,
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
        "max_frobenius_error": s.max_frobenius_error,
        "mean_frobenius_error": s.mean_frobenius_error,
        "max_rotation_angle_rad": s.max_rotation_angle,
        "mean_rotation_angle_rad": s.mean_rotation_angle,
        "max_unitarity_error": s.max_unitarity_error,
        "max_det_error": s.max_det_error,
    }


def main() -> None:
    n = 10
    _, U = random_pure_gauge(n, seed=101)
    H = random_su2_stack(n, seed=102)
    U_gauge = gauge_transform(U, H)

    corrupt_edges = [(0, 1), (2, 5), (3, 7), (6, 9)]
    U_bad, metadata = perturb_links(U, corrupt_edges, strength=1.0, seed=103)
    U_bad_gauge = gauge_transform(U_bad, H)

    rows = [
        summary_row("pure_gauge", U),
        summary_row("pure_gauge_after_local_gauge_transform", U_gauge),
        summary_row("perturbed_sparse_links", U_bad),
        summary_row("perturbed_after_local_gauge_transform", U_bad_gauge),
    ]

    cov_err_pure = residual_covariance_error(U, H)
    cov_err_bad = residual_covariance_error(U_bad, H)
    literal_change_bad = residual_literal_change(U_bad, H)
    trace_inv_bad = residual_trace_invariance_error(U_bad, H)

    csv_path = OUT_TABLES / "su2_nonabelian_covariance_summary.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
        writer.writerow({})
        writer.writerow({"case": "covariance_error_pure", "max_frobenius_error": cov_err_pure})
        writer.writerow({"case": "covariance_error_perturbed", "max_frobenius_error": cov_err_bad})
        writer.writerow({"case": "literal_residual_matrix_change_perturbed", "max_frobenius_error": literal_change_bad})
        writer.writerow({"case": "trace_invariance_error_perturbed", "max_frobenius_error": trace_inv_bad})
        for edge, meta in metadata.items():
            writer.writerow({"case": f"corruption_edge_{edge}", "max_rotation_angle_rad": meta["angle"]})

    _, R_bad = triangle_residuals(U_bad)
    angles = su2_rotation_angles(R_bad)
    fig_path = OUT_FIGS / "su2_perturbed_residual_histogram.png"
    plt.figure(figsize=(8, 5))
    plt.hist(angles, bins=40)
    plt.xlabel("SU(2) residual rotation angle |theta| [rad]")
    plt.ylabel("ordered triangle count")
    plt.title("SU(2) Her re-anchoring residuals after sparse link perturbations")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=160)
    plt.close()

    scores = residual_participation_scores(U_bad)
    top_m = top_scored_edges(scores, len(corrupt_edges))
    top_2m = top_scored_edges(scores, 2 * len(corrupt_edges))
    true_edges = {tuple(sorted(e)) for e in corrupt_edges}
    recall_m = len({e for e, _ in top_m} & true_edges)
    recall_2m = len({e for e, _ in top_2m} & true_edges)

    print("SU(2) non-Abelian covariance experiment")
    print(f"pure gauge max residual:       {rows[0]['max_frobenius_error']:.3e}")
    print(f"perturbed max residual:        {rows[2]['max_frobenius_error']:.3e}")
    print(f"covariance error, pure:        {cov_err_pure:.3e}")
    print(f"covariance error, perturbed:   {cov_err_bad:.3e}")
    print(f"literal residual matrix change:{literal_change_bad:.3e}")
    print(f"trace invariance error:        {trace_inv_bad:.3e}")
    print(f"top-m recall:                  {recall_m}/{len(corrupt_edges)}")
    print(f"top-2m recall:                 {recall_2m}/{len(corrupt_edges)}")
    print(f"wrote {csv_path}")
    print(f"wrote {fig_path}")


if __name__ == "__main__":
    main()
