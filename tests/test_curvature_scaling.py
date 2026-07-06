import numpy as np

from her_reanchoring.curvature_scaling import (
    constant_curvature_field,
    constant_curvature_potential,
    curvature_scaling_table,
    curvature_values_at_centers,
    estimate_loglog_slope,
    her_reconstructed_angle_density,
    max_u1_plaquette_gauge_invariance_error_for_smooth_field,
    plaquette_angle_density,
    u1_midpoint_link_matrix,
)


def test_constant_u1_curvature_scales_exactly_with_area():
    B = 0.37
    potential = constant_curvature_potential(B)
    curvature = constant_curvature_field(B)
    for n_side in (4, 7, 11):
        U = u1_midpoint_link_matrix(n_side, potential)
        density = plaquette_angle_density(U, n_side)
        target = curvature_values_at_centers(n_side, curvature)
        assert np.max(np.abs(density - target)) < 1e-12


def test_her_reconstructed_density_matches_ordinary_plaquette_density():
    n_side = 9
    U = u1_midpoint_link_matrix(n_side, constant_curvature_potential(0.23))
    ordinary = plaquette_angle_density(U, n_side)
    her = her_reconstructed_angle_density(U, n_side)
    assert np.max(np.abs(ordinary - her)) < 1e-12


def test_smooth_u1_curvature_error_decreases_under_refinement():
    rows = curvature_scaling_table((8, 12, 16, 24, 32))
    errors = [row.rms_curvature_error for row in rows]
    hs = [row.h for row in rows]
    assert errors[-1] < 0.2 * errors[0]
    assert estimate_loglog_slope(hs, errors) > 1.7
    assert max(row.max_her_reconstruction_error for row in rows) < 1e-12


def test_smooth_u1_plaquettes_remain_gauge_invariant():
    n_side = 8
    rng = np.random.default_rng(123)
    lambdas = rng.uniform(-np.pi, np.pi, size=n_side * n_side)
    err = max_u1_plaquette_gauge_invariance_error_for_smooth_field(n_side, lambdas)
    assert err < 1e-12
