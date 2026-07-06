import numpy as np

from her_reanchoring.su2_gauge import (
    I2,
    dagger,
    det_errors,
    gauge_transform,
    perturb_links,
    random_pure_gauge,
    random_su2_stack,
    reanchored_link_covariance_error,
    residual_covariance_error,
    residual_literal_change,
    residual_participation_scores,
    residual_summary,
    residual_trace_invariance_error,
    top_scored_edges,
    triangle_residual,
    unitarity_errors,
)


def test_random_su2_frames_are_unitary_with_determinant_one():
    H = random_su2_stack(12, seed=1)
    assert np.max(unitarity_errors(H)) < 1e-12
    assert np.max(det_errors(H)) < 1e-12


def test_su2_pure_gauge_triangle_residuals_are_machine_zero():
    _, U = random_pure_gauge(7, seed=2)
    summary = residual_summary(U)
    assert summary.max_frobenius_error < 1e-12
    assert summary.max_rotation_angle < 1e-12
    assert summary.max_unitarity_error < 1e-12
    assert summary.max_det_error < 1e-12


def test_su2_gauge_transform_preserves_pure_gauge_residual_zero():
    _, U = random_pure_gauge(7, seed=3)
    H = random_su2_stack(7, seed=4)
    V = gauge_transform(U, H)
    summary = residual_summary(V)
    assert summary.max_frobenius_error < 1e-12
    assert summary.max_rotation_angle < 1e-12


def test_su2_reanchored_prediction_transforms_as_j_to_k_link():
    _, U = random_pure_gauge(6, seed=5)
    H = random_su2_stack(6, seed=6)
    err = reanchored_link_covariance_error(U, H)
    assert err < 1e-12


def test_su2_residuals_are_covariant_not_literally_invariant():
    _, U = random_pure_gauge(8, seed=7)
    V, _ = perturb_links(U, [(0, 1), (2, 5), (3, 7)], strength=1.1, seed=8)
    H = random_su2_stack(8, seed=9)

    cov_err = residual_covariance_error(V, H)
    literal_change = residual_literal_change(V, H)
    trace_err = residual_trace_invariance_error(V, H)

    assert cov_err < 1e-12
    assert trace_err < 1e-12
    assert literal_change > 1e-3


def test_perturbed_su2_links_create_nonzero_residuals():
    _, U = random_pure_gauge(8, seed=10)
    V, _ = perturb_links(U, [(0, 1)], strength=0.8, seed=11)
    summary = residual_summary(V)
    assert summary.max_frobenius_error > 1e-3
    assert summary.mean_rotation_angle > 1e-4


def test_sparse_fault_scores_rank_true_edges_highly_su2():
    _, U = random_pure_gauge(10, seed=12)
    true_edges = {(0, 1), (2, 5), (3, 7)}
    V, _ = perturb_links(U, list(true_edges), strength=0.9, seed=13)
    scores = residual_participation_scores(V)
    top = {edge for edge, _ in top_scored_edges(scores, len(true_edges))}
    assert len(top & true_edges) >= 2


def test_triangle_residual_identity_for_single_triangle():
    _, U = random_pure_gauge(5, seed=14)
    R = triangle_residual(U, 0, 1, 2)
    assert np.linalg.norm(R - I2) < 1e-12
    assert np.linalg.norm(dagger(R) @ R - I2) < 1e-12
