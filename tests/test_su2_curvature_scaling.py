import numpy as np

from her_reanchoring.su2_curvature_scaling import (
    abelian_part_of_default_su2_curvature,
    default_su2_curvature_vector,
    estimate_loglog_slope,
    max_su2_her_reconstruction_error_for_smooth_field,
    max_su2_plaquette_covariance_error_for_smooth_field,
    max_su2_plaquette_trace_invariance_error_for_smooth_field,
    su2_curvature_scaling_table,
    su2_exp_from_vector,
    su2_her_reconstructed_log_density,
    su2_log_vector,
    su2_midpoint_link_matrix,
    su2_plaquette_log_density,
)


def test_su2_exp_log_round_trip_near_identity():
    vectors = [
        np.array([0.01, -0.02, 0.03]),
        np.array([0.12, 0.04, -0.07]),
        np.array([-0.09, 0.08, 0.02]),
    ]
    for v in vectors:
        assert np.linalg.norm(su2_log_vector(su2_exp_from_vector(v)) - v) < 1e-12


def test_default_su2_curvature_has_nonzero_commutator_term():
    x = np.array([0.25, 0.75])
    y = np.array([0.30, 0.60])
    full = default_su2_curvature_vector(x, y)
    abelian_only = abelian_part_of_default_su2_curvature(x, y)
    assert np.max(np.linalg.norm(full - abelian_only, axis=1)) > 1e-3


def test_her_reconstructed_su2_density_matches_ordinary_plaquette_density():
    n_side = 8
    U = su2_midpoint_link_matrix(n_side)
    ordinary = su2_plaquette_log_density(U, n_side)
    her = su2_her_reconstructed_log_density(U, n_side)
    assert np.max(np.linalg.norm(ordinary - her, axis=1)) < 1e-11


def test_smooth_su2_curvature_error_decreases_under_refinement():
    rows = su2_curvature_scaling_table((6, 8, 10, 12, 16, 20))
    vector_errors = [row.rms_vector_curvature_error for row in rows]
    norm_errors = [row.rms_norm_curvature_error for row in rows]
    hs = [row.h for row in rows]

    # In a non-Abelian field, log(W) is based at the plaquette basepoint while
    # the analytic F target is sampled at the center. Without explicit parallel
    # transport, the raw vector comparison is expected to converge, but only
    # roughly first order. Conjugacy-invariant norm diagnostics recover the
    # second-order midpoint scaling.
    assert vector_errors[-1] < 0.35 * vector_errors[0]
    assert 0.7 < estimate_loglog_slope(hs, vector_errors) < 1.3
    assert norm_errors[-1] < 0.08 * norm_errors[0]
    assert estimate_loglog_slope(hs, norm_errors) > 1.7
    assert max(row.max_her_reconstruction_error for row in rows) < 1e-11


def test_smooth_su2_plaquettes_transform_covariantly():
    assert max_su2_plaquette_covariance_error_for_smooth_field(7, seed=321) < 1e-12
    assert max_su2_plaquette_trace_invariance_error_for_smooth_field(7, seed=321) < 1e-12


def test_smooth_su2_her_reconstruction_is_machine_precision():
    assert max_su2_her_reconstruction_error_for_smooth_field(9) < 1e-12
