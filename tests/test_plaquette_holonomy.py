import numpy as np

from her_reanchoring.plaquette_holonomy import (
    her_plaquette_reconstruction_error,
    her_reconstructed_plaquette_holonomy,
    her_triangle_holonomy_conjugacy_error,
    plaquette_holonomy,
    plaquette_summary,
    square_plaquettes,
    su2_plaquette_covariance_error,
    su2_plaquette_trace_invariance_error,
    triangle_factorization_error,
    u1_plaquette_gauge_invariance_error,
)
from her_reanchoring.u1_gauge import perturb_links as perturb_u1_links
from her_reanchoring.u1_gauge import random_pure_gauge as random_u1_pure_gauge
from her_reanchoring.su2_gauge import (
    random_pure_gauge as random_su2_pure_gauge,
    random_su2_stack,
    perturb_links as perturb_su2_links,
)


def test_u1_pure_gauge_plaquettes_are_flat():
    _, U = random_u1_pure_gauge(9, seed=10)
    summary = plaquette_summary(U, nx=3, ny=3, group="u1")
    assert summary.n_plaquettes == 4
    assert summary.max_plaquette_angle < 1e-12
    assert summary.max_her_reconstruction_error < 1e-12


def test_u1_plaquette_equals_two_triangles_and_her_reconstruction():
    _, U0 = random_u1_pure_gauge(9, seed=11)
    U, _ = perturb_u1_links(U0, edges=[(1, 2), (4, 7)], strength=0.7, seed=12)
    for p in square_plaquettes(3, 3):
        assert triangle_factorization_error(U, p, "u1") < 1e-12
        assert her_plaquette_reconstruction_error(U, p, "u1") < 1e-12


def test_u1_plaquettes_are_gauge_invariant():
    _, U0 = random_u1_pure_gauge(9, seed=13)
    U, _ = perturb_u1_links(U0, edges=[(0, 1), (5, 8)], strength=0.8, seed=14)
    rng = np.random.default_rng(15)
    lambdas = rng.uniform(-np.pi, np.pi, size=9)
    assert u1_plaquette_gauge_invariance_error(U, lambdas, nx=3, ny=3) < 1e-12


def test_su2_pure_gauge_plaquettes_are_flat():
    _, U = random_su2_pure_gauge(9, seed=20)
    summary = plaquette_summary(U, nx=3, ny=3, group="su2")
    assert summary.n_plaquettes == 4
    assert summary.max_plaquette_angle < 1e-12
    assert summary.max_her_reconstruction_error < 1e-12


def test_su2_plaquette_equals_two_triangles_and_her_reconstruction():
    _, U0 = random_su2_pure_gauge(9, seed=21)
    U, _ = perturb_su2_links(U0, edges=[(1, 2), (4, 7)], strength=0.7, seed=22)
    for p in square_plaquettes(3, 3):
        assert triangle_factorization_error(U, p, "su2") < 1e-12
        assert her_plaquette_reconstruction_error(U, p, "su2") < 1e-12


def test_su2_her_triangle_residual_is_conjugate_to_inverse_holonomy():
    _, U0 = random_su2_pure_gauge(9, seed=23)
    U, _ = perturb_su2_links(U0, edges=[(1, 2), (4, 7)], strength=0.7, seed=24)
    for p in square_plaquettes(3, 3):
        a, b, c, d = p
        assert her_triangle_holonomy_conjugacy_error(U, a, b, c, "su2") < 1e-12
        assert her_triangle_holonomy_conjugacy_error(U, a, c, d, "su2") < 1e-12


def test_su2_plaquettes_transform_covariantly_not_literally():
    _, U0 = random_su2_pure_gauge(9, seed=25)
    U, _ = perturb_su2_links(U0, edges=[(1, 2), (4, 7)], strength=0.9, seed=26)
    H = random_su2_stack(9, seed=27)
    assert su2_plaquette_covariance_error(U, H, nx=3, ny=3) < 1e-12
    assert su2_plaquette_trace_invariance_error(U, H, nx=3, ny=3) < 1e-12


def test_her_reconstruction_matches_plaquette_even_when_not_flat():
    _, U0 = random_su2_pure_gauge(9, seed=30)
    U, _ = perturb_su2_links(U0, edges=[(0, 1), (1, 2), (3, 6)], strength=1.0, seed=31)
    differences = []
    nonflat_angles = []
    for p in square_plaquettes(3, 3):
        W = plaquette_holonomy(U, p, "su2")
        W_her = her_reconstructed_plaquette_holonomy(U, p, "su2")
        differences.append(np.linalg.norm(W - W_her))
        nonflat_angles.append(np.linalg.norm(W - np.eye(2)))
    assert max(differences) < 1e-12
    assert max(nonflat_angles) > 1e-3
