"""Core Her Algorithm physics helpers."""

from .u1_gauge import (
    ResidualSummary,
    angle_errors,
    as_u1,
    chord_errors,
    gauge_residual_covariance_error,
    gauge_transform,
    inverse_u1,
    perturb_links,
    pure_gauge_links,
    random_pure_gauge,
    random_reciprocal_links,
    reanchored_link,
    residual_participation_scores,
    residual_summary,
    top_scored_edges,
    triangle_residual,
    triangle_residuals,
)

__all__ = [
    "ResidualSummary",
    "angle_errors",
    "as_u1",
    "chord_errors",
    "gauge_residual_covariance_error",
    "gauge_transform",
    "inverse_u1",
    "perturb_links",
    "pure_gauge_links",
    "random_pure_gauge",
    "random_reciprocal_links",
    "reanchored_link",
    "residual_participation_scores",
    "residual_summary",
    "top_scored_edges",
    "triangle_residual",
    "triangle_residuals",
]


from .su2_gauge import (
    MatrixResidualSummary,
    dagger,
    det_errors,
    random_su2,
    random_su2_stack,
    residual_covariance_error as su2_residual_covariance_error,
    residual_literal_change as su2_residual_literal_change,
    residual_trace_invariance_error as su2_residual_trace_invariance_error,
    su2_axis_angle,
    su2_rotation_angles,
    unitarity_errors,
)

__all__ += [
    "MatrixResidualSummary",
    "dagger",
    "det_errors",
    "random_su2",
    "random_su2_stack",
    "su2_residual_covariance_error",
    "su2_residual_literal_change",
    "su2_residual_trace_invariance_error",
    "su2_axis_angle",
    "su2_rotation_angles",
    "unitarity_errors",
]
