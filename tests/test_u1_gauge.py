import numpy as np

from her_reanchoring.u1_gauge import (
    gauge_residual_covariance_error,
    gauge_transform,
    perturb_links,
    pure_gauge_links,
    random_pure_gauge,
    residual_participation_scores,
    residual_summary,
    top_scored_edges,
    triangle_residuals,
)


def test_pure_gauge_triangle_residuals_are_machine_zero():
    _, U = random_pure_gauge(8, seed=1)
    summary = residual_summary(U)
    assert summary.max_chord_error < 1e-12
    assert summary.max_angle_error < 1e-12


def test_u1_gauge_transform_preserves_pure_gauge_residuals():
    _, U = random_pure_gauge(8, seed=2)
    lambdas = np.linspace(0.1, 1.3, 8)
    V = gauge_transform(U, lambdas)
    summary = residual_summary(V)
    assert summary.max_chord_error < 1e-12
    assert summary.max_angle_error < 1e-12


def test_u1_residuals_are_gauge_invariant_for_perturbed_field():
    _, U = random_pure_gauge(9, seed=3)
    V, _ = perturb_links(U, [(0, 1), (2, 5), (3, 7)], strength=0.8, seed=4)
    lambdas = np.random.default_rng(5).uniform(-np.pi, np.pi, size=9)
    err = gauge_residual_covariance_error(V, lambdas)
    assert err < 1e-12


def test_perturbed_links_create_nonzero_residuals():
    _, U = random_pure_gauge(8, seed=6)
    V, _ = perturb_links(U, [(0, 1)], strength=0.5, seed=7)
    summary = residual_summary(V)
    assert summary.max_angle_error > 1e-3
    assert summary.mean_angle_error > 1e-4


def test_sparse_fault_scores_rank_true_edges_highly():
    _, U = random_pure_gauge(10, seed=8)
    true_edges = {(0, 1), (2, 5), (3, 7)}
    V, _ = perturb_links(U, list(true_edges), strength=0.8, seed=9)
    scores = residual_participation_scores(V)
    top = {edge for edge, _ in top_scored_edges(scores, len(true_edges))}
    assert len(top & true_edges) >= 2
