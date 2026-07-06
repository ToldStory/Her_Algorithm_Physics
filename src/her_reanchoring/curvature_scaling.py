"""Curvature scaling tests for U(1) smooth lattice connections.

This module is deliberately conservative.  Earlier experiments proved exact
algebraic compatibility between Her triangle residuals and plaquette holonomy.
Here we ask a different question:

    If lattice links are sampled from a smooth U(1) gauge potential A, does the
    Her-reconstructed plaquette holonomy approach exp(i F_xy h^2) as the lattice
    spacing h shrinks?

For U(1), the continuum curvature is

    F_xy = partial_x A_y - partial_y A_x.

The midpoint links below approximate line integrals along straight segments.
For a square plaquette of side h, the plaquette angle divided by h^2 should
converge to F_xy at the plaquette center.  The Her-reconstructed plaquette is
then tested against the same continuum target through the exact reconstruction
identity already established in the plaquette/holonomy layer.

This is not a claim that Her residuals are identical to curvature.  It is a
scaling-limit test: Her residuals reconstruct plaquette holonomy, and plaquette
holonomy has the expected small-area curvature behavior in this synthetic U(1)
case.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

import numpy as np

from .plaquette_holonomy import (
    her_plaquette_reconstruction_error,
    her_reconstructed_plaquette_holonomy,
    plaquette_holonomies,
    square_plaquettes,
)
from .u1_gauge import as_u1, gauge_transform as u1_gauge_transform

Array = np.ndarray
ScalarField = Callable[[Array | float, Array | float], Array | float]
VectorPotential = Callable[[Array | float, Array | float], tuple[Array | float, Array | float]]


@dataclass(frozen=True)
class CurvatureScalingRow:
    """One grid-resolution row for the curvature scaling experiment."""

    n_side: int
    h: float
    n_plaquettes: int
    rms_curvature_error: float
    max_abs_curvature_error: float
    rms_angle_error: float
    max_abs_angle_error: float
    max_her_reconstruction_error: float


def default_u1_vector_potential(
    x: Array | float,
    y: Array | float,
) -> tuple[Array | float, Array | float]:
    """A smooth synthetic U(1) gauge potential on the unit square.

    The coefficients are intentionally modest so plaquette angles remain far
    from the branch cut of the wrapped phase during refinement tests.
    """

    two_pi = 2.0 * np.pi
    ax = 0.19 * np.sin(two_pi * y) + 0.07 * x * y
    ay = -0.13 * np.cos(two_pi * x) + 0.05 * x * x * y + 0.11 * x
    return ax, ay


def default_u1_curvature(x: Array | float, y: Array | float) -> Array | float:
    """Analytic curvature for :func:`default_u1_vector_potential`.

    F_xy = partial_x A_y - partial_y A_x.
    """

    two_pi = 2.0 * np.pi
    partial_x_ay = 0.13 * two_pi * np.sin(two_pi * x) + 0.10 * x * y + 0.11
    partial_y_ax = 0.19 * two_pi * np.cos(two_pi * y) + 0.07 * x
    return partial_x_ay - partial_y_ax


def constant_curvature_potential(B: float) -> VectorPotential:
    """Return A=(0, B x), whose continuum curvature is exactly B."""

    def potential(x: Array | float, y: Array | float) -> tuple[Array | float, Array | float]:
        return np.zeros_like(np.asarray(x, dtype=float)), B * np.asarray(x, dtype=float)

    return potential


def constant_curvature_field(B: float) -> ScalarField:
    """Return the constant scalar curvature F_xy=B."""

    def curvature(x: Array | float, y: Array | float) -> Array | float:
        return np.zeros_like(np.asarray(x, dtype=float)) + B

    return curvature


def grid_positions(n_side: int) -> Array:
    """Return row-major (x,y) site positions on [0,1]^2."""

    if n_side < 2:
        raise ValueError("n_side must be at least 2.")
    xs = np.linspace(0.0, 1.0, n_side)
    ys = np.linspace(0.0, 1.0, n_side)
    return np.array([(x, y) for y in ys for x in xs], dtype=float)


def u1_midpoint_link_matrix(n_side: int, potential: VectorPotential) -> Array:
    """Build complete directed U(1) links from midpoint line-integral samples.

    For sites p_i and p_j, the link is

        U_ij = exp(i (p_j-p_i) dot A((p_i+p_j)/2)).

    This is exact for constant vector potentials along a segment and a
    second-order local approximation for smooth potentials.  Complete links are
    generated because Her re-anchoring uses triangle diagonals in addition to
    nearest-neighbor plaquette edges.
    """

    positions = grid_positions(n_side)
    n_sites = positions.shape[0]
    U = np.ones((n_sites, n_sites), dtype=complex)
    for i in range(n_sites):
        pi = positions[i]
        for j in range(i + 1, n_sites):
            pj = positions[j]
            delta = pj - pi
            mid = 0.5 * (pi + pj)
            ax, ay = potential(mid[0], mid[1])
            angle = float(delta[0] * ax + delta[1] * ay)
            U[i, j] = as_u1(angle)
            U[j, i] = np.conjugate(U[i, j])
    return U


def plaquette_centers(n_side: int) -> Array:
    """Return row-major square plaquette centers on [0,1]^2."""

    if n_side < 2:
        raise ValueError("n_side must be at least 2.")
    h = 1.0 / (n_side - 1)
    centers = []
    for y in range(n_side - 1):
        for x in range(n_side - 1):
            centers.append(((x + 0.5) * h, (y + 0.5) * h))
    return np.array(centers, dtype=float)


def wrapped_angle(z: Array) -> Array:
    """Return principal U(1) phase angles in (-pi, pi]."""

    return np.angle(z)


def plaquette_angle_density(U: Array, n_side: int) -> Array:
    """Return plaquette angles divided by h^2."""

    h = 1.0 / (n_side - 1)
    plaquettes = square_plaquettes(n_side, n_side)
    W = plaquette_holonomies(U, plaquettes, "u1")
    return wrapped_angle(W) / (h * h)


def her_reconstructed_angle_density(U: Array, n_side: int) -> Array:
    """Return Her-reconstructed plaquette angles divided by h^2."""

    h = 1.0 / (n_side - 1)
    plaquettes = square_plaquettes(n_side, n_side)
    W_her = np.array([her_reconstructed_plaquette_holonomy(U, p, "u1") for p in plaquettes])
    return wrapped_angle(W_her) / (h * h)


def curvature_values_at_centers(n_side: int, curvature: ScalarField) -> Array:
    """Evaluate analytic curvature at square plaquette centers."""

    centers = plaquette_centers(n_side)
    return np.asarray(curvature(centers[:, 0], centers[:, 1]), dtype=float)


def curvature_scaling_row(
    n_side: int,
    potential: VectorPotential = default_u1_vector_potential,
    curvature: ScalarField = default_u1_curvature,
) -> CurvatureScalingRow:
    """Compute one curvature-scaling summary row."""

    U = u1_midpoint_link_matrix(n_side, potential)
    h = 1.0 / (n_side - 1)
    density = plaquette_angle_density(U, n_side)
    her_density = her_reconstructed_angle_density(U, n_side)
    target = curvature_values_at_centers(n_side, curvature)
    err = density - target
    angle_err = wrapped_angle(np.exp(1j * (density - target) * h * h))
    plaquettes = square_plaquettes(n_side, n_side)
    her_errors = np.array([her_plaquette_reconstruction_error(U, p, "u1") for p in plaquettes])
    her_density_err = her_density - density
    return CurvatureScalingRow(
        n_side=n_side,
        h=float(h),
        n_plaquettes=len(plaquettes),
        rms_curvature_error=float(np.sqrt(np.mean(err**2))),
        max_abs_curvature_error=float(np.max(np.abs(err))),
        rms_angle_error=float(np.sqrt(np.mean(angle_err**2))),
        max_abs_angle_error=float(np.max(np.abs(angle_err))),
        max_her_reconstruction_error=float(max(np.max(her_errors), np.max(np.abs(her_density_err)) * h * h)),
    )


def curvature_scaling_table(
    n_sides: Iterable[int] = (8, 12, 16, 24, 32),
    potential: VectorPotential = default_u1_vector_potential,
    curvature: ScalarField = default_u1_curvature,
) -> list[CurvatureScalingRow]:
    """Compute scaling rows for several grid resolutions."""

    return [curvature_scaling_row(int(n), potential, curvature) for n in n_sides]


def estimate_loglog_slope(h_values: Sequence[float], errors: Sequence[float]) -> float:
    """Estimate slope of log(error) versus log(h)."""

    h = np.asarray(h_values, dtype=float)
    e = np.asarray(errors, dtype=float)
    mask = (h > 0) & (e > 0)
    if int(np.sum(mask)) < 2:
        raise ValueError("Need at least two positive h/error pairs to estimate a slope.")
    slope, _ = np.polyfit(np.log(h[mask]), np.log(e[mask]), deg=1)
    return float(slope)


def max_u1_plaquette_gauge_invariance_error_for_smooth_field(
    n_side: int,
    lambdas: Sequence[float],
    potential: VectorPotential = default_u1_vector_potential,
) -> float:
    """Check that plaquette holonomies remain invariant under local U(1) gauge changes."""

    U = u1_midpoint_link_matrix(n_side, potential)
    V = u1_gauge_transform(U, lambdas)
    plaquettes = square_plaquettes(n_side, n_side)
    W = plaquette_holonomies(U, plaquettes, "u1")
    W2 = plaquette_holonomies(V, plaquettes, "u1")
    return float(np.max(np.abs(W2 - W))) if len(plaquettes) else 0.0
