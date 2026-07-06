"""Non-Abelian SU(2) curvature scaling tests.

This module extends the U(1) curvature-scaling experiment to a smooth
synthetic SU(2) connection.  It is deliberately conservative: the tested claim
is not that Her residuals *are* curvature.  The tested claim is that Her
residuals reconstruct plaquette holonomy, and that the reconstructed plaquette
holonomy has the expected small-area relationship to the non-Abelian curvature
of a smooth SU(2) connection.

Conventions
-----------
We represent a Hermitian SU(2) gauge potential as

    A_x(x,y) = a_x(x,y) . sigma
    A_y(x,y) = a_y(x,y) . sigma,

where ``sigma`` are the Pauli matrices and ``a_x, a_y`` are real 3-vectors.
Links are sampled by a midpoint approximation

    U_ij = exp(i ((p_j-p_i)_x A_x(mid) + (p_j-p_i)_y A_y(mid))).

With this convention, the Hermitian non-Abelian curvature vector is

    F_xy = partial_x a_y - partial_y a_x - 2 a_x cross a_y,

because ``i[A_x,A_y] = -2 (a_x cross a_y).sigma``.

For a small square of area h^2,

    W_square ~= exp(i F_xy h^2),

so the principal SU(2) logarithm vector divided by h^2 should converge toward
the analytic curvature vector.  A subtle non-Abelian basepoint issue remains:
plaquette logs are based at a plaquette corner, while the analytic target is
sampled at the center.  Without explicitly parallel-transporting the target
into the same fiber, the raw vector comparison converges roughly first order;
conjugacy-invariant norm diagnostics recover the expected second-order midpoint
scaling.  The comparison is made only in a modest-amplitude regime, away from
the SU(2) logarithm branch cut.
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
    su2_plaquette_covariance_error,
    su2_plaquette_trace_invariance_error,
)
from .su2_gauge import I2, PAULI, gauge_transform as su2_gauge_transform, random_su2_stack

Array = np.ndarray
Vector3Field = Callable[[Array | float, Array | float], Array]
SU2Potential = Callable[[Array | float, Array | float], tuple[Array, Array]]


@dataclass(frozen=True)
class SU2CurvatureScalingRow:
    """One grid-resolution row for the SU(2) curvature scaling experiment."""

    n_side: int
    h: float
    n_plaquettes: int
    rms_vector_curvature_error: float
    max_vector_curvature_error: float
    rms_norm_curvature_error: float
    max_norm_curvature_error: float
    max_her_reconstruction_error: float
    max_plaquette_log_norm: float


def _as_broadcast_arrays(x: Array | float, y: Array | float) -> tuple[Array, Array]:
    """Return x,y as broadcast-compatible float arrays."""

    return np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))


def _stack3(a: Array, b: Array, c: Array) -> Array:
    """Stack three scalar arrays into a final-axis vector field."""

    return np.stack(np.broadcast_arrays(a, b, c), axis=-1).astype(float)


def su2_exp_from_vector(v: Sequence[float] | Array) -> Array:
    """Return exp(i v.sigma) as an SU(2) matrix.

    ``v`` is a real 3-vector.  For small vectors this is ``I + i v.sigma`` to
    first order.  The exact closed form is used to avoid scipy dependencies.
    """

    vec = np.asarray(v, dtype=float)
    if vec.shape != (3,):
        raise ValueError("su2_exp_from_vector expects a length-3 vector.")
    r = float(np.linalg.norm(vec))
    if r < 1e-14:
        return I2 + 1j * (vec[0] * PAULI[0] + vec[1] * PAULI[1] + vec[2] * PAULI[2])
    n = vec / r
    generator = n[0] * PAULI[0] + n[1] * PAULI[1] + n[2] * PAULI[2]
    return np.cos(r) * I2 + 1j * np.sin(r) * generator


def su2_log_vector(U: Array) -> Array:
    """Principal logarithm vector v for an SU(2) matrix U ~= exp(i v.sigma).

    The return vector has norm in [0, pi].  The refinement tests use matrices
    close to identity, so branch ambiguity is intentionally avoided.
    """

    M = np.asarray(U, dtype=complex)
    if M.shape != (2, 2):
        raise ValueError("su2_log_vector expects one 2x2 matrix.")
    cos_r = float(np.clip(np.real(np.trace(M)) / 2.0, -1.0, 1.0))
    r = float(np.arccos(cos_r))
    coeff = np.array(
        [np.real(-0.5j * np.trace(PAULI[a] @ M)) for a in range(3)],
        dtype=float,
    )
    sin_r = float(np.sin(r))
    if abs(sin_r) < 1e-12:
        # Near the identity, coeff = sin(r) n ~= r n = v.
        return coeff
    return (r / sin_r) * coeff


def su2_log_vectors(mats: Array) -> Array:
    """Apply :func:`su2_log_vector` to a stack of SU(2) matrices."""

    M = np.asarray(mats, dtype=complex)
    return np.array([su2_log_vector(U) for U in M.reshape((-1, 2, 2))]).reshape(M.shape[:-2] + (3,))


def default_su2_vector_potential(
    x: Array | float,
    y: Array | float,
) -> tuple[Array, Array]:
    """A smooth non-Abelian synthetic SU(2) potential.

    The amplitude is intentionally small enough that plaquette logarithms stay
    away from branch cuts, while still including a nonzero commutator term.
    """

    x_arr, y_arr = _as_broadcast_arrays(x, y)
    ax = _stack3(
        0.00 * x_arr,
        0.17 * y_arr,
        0.03 * np.sin(2.0 * np.pi * x_arr),
    )
    ay = _stack3(
        0.13 * x_arr,
        0.02 * np.cos(2.0 * np.pi * y_arr),
        0.04 * x_arr * y_arr,
    )
    return ax, ay


def default_su2_curvature_vector(x: Array | float, y: Array | float) -> Array:
    """Analytic non-Abelian curvature vector for default potential.

    F_xy = partial_x a_y - partial_y a_x - 2 a_x cross a_y.
    """

    x_arr, y_arr = _as_broadcast_arrays(x, y)
    ax, ay = default_su2_vector_potential(x_arr, y_arr)
    partial_x_ay = _stack3(
        np.zeros_like(x_arr) + 0.13,
        np.zeros_like(x_arr),
        0.04 * y_arr,
    )
    partial_y_ax = _stack3(
        np.zeros_like(x_arr),
        np.zeros_like(x_arr) + 0.17,
        np.zeros_like(x_arr),
    )
    return partial_x_ay - partial_y_ax - 2.0 * np.cross(ax, ay)


def abelian_part_of_default_su2_curvature(x: Array | float, y: Array | float) -> Array:
    """Derivative-only part of the default SU(2) curvature.

    This is included only for tests/documentation showing that the default case
    is genuinely non-Abelian; the true target includes the commutator term.
    """

    x_arr, y_arr = _as_broadcast_arrays(x, y)
    partial_x_ay = _stack3(
        np.zeros_like(x_arr) + 0.13,
        np.zeros_like(x_arr),
        0.04 * y_arr,
    )
    partial_y_ax = _stack3(
        np.zeros_like(x_arr),
        np.zeros_like(x_arr) + 0.17,
        np.zeros_like(x_arr),
    )
    return partial_x_ay - partial_y_ax


def su2_grid_positions(n_side: int) -> Array:
    """Return row-major (x,y) site positions on [0,1]^2."""

    if n_side < 2:
        raise ValueError("n_side must be at least 2.")
    xs = np.linspace(0.0, 1.0, n_side)
    ys = np.linspace(0.0, 1.0, n_side)
    return np.array([(x, y) for y in ys for x in xs], dtype=float)


def su2_midpoint_link_matrix(
    n_side: int,
    potential: SU2Potential = default_su2_vector_potential,
) -> Array:
    """Build complete directed SU(2) links by midpoint line-integral sampling."""

    positions = su2_grid_positions(n_side)
    n_sites = positions.shape[0]
    U = np.empty((n_sites, n_sites, 2, 2), dtype=complex)
    for i in range(n_sites):
        U[i, i] = I2
    for i in range(n_sites):
        pi = positions[i]
        for j in range(i + 1, n_sites):
            pj = positions[j]
            delta = pj - pi
            mid = 0.5 * (pi + pj)
            ax, ay = potential(mid[0], mid[1])
            # Squeeze handles scalar-input vector fields returning shape (3,) or (1,3).
            vec = delta[0] * np.squeeze(ax) + delta[1] * np.squeeze(ay)
            U[i, j] = su2_exp_from_vector(vec)
            U[j, i] = U[i, j].conj().T
    return U


def su2_plaquette_centers(n_side: int) -> Array:
    """Return row-major square plaquette centers on [0,1]^2."""

    if n_side < 2:
        raise ValueError("n_side must be at least 2.")
    h = 1.0 / (n_side - 1)
    centers = []
    for y in range(n_side - 1):
        for x in range(n_side - 1):
            centers.append(((x + 0.5) * h, (y + 0.5) * h))
    return np.array(centers, dtype=float)


def su2_curvature_values_at_centers(
    n_side: int,
    curvature: Vector3Field = default_su2_curvature_vector,
) -> Array:
    """Evaluate analytic SU(2) curvature vectors at plaquette centers."""

    centers = su2_plaquette_centers(n_side)
    return np.asarray(curvature(centers[:, 0], centers[:, 1]), dtype=float)


def su2_plaquette_log_density(U: Array, n_side: int) -> Array:
    """Return principal plaquette log vectors divided by h^2."""

    h = 1.0 / (n_side - 1)
    plaquettes = square_plaquettes(n_side, n_side)
    W = plaquette_holonomies(U, plaquettes, "su2")
    return su2_log_vectors(W) / (h * h)


def su2_her_reconstructed_log_density(U: Array, n_side: int) -> Array:
    """Return Her-reconstructed plaquette log vectors divided by h^2."""

    h = 1.0 / (n_side - 1)
    plaquettes = square_plaquettes(n_side, n_side)
    W_her = np.array([her_reconstructed_plaquette_holonomy(U, p, "su2") for p in plaquettes])
    return su2_log_vectors(W_her) / (h * h)


def su2_curvature_scaling_row(
    n_side: int,
    potential: SU2Potential = default_su2_vector_potential,
    curvature: Vector3Field = default_su2_curvature_vector,
) -> SU2CurvatureScalingRow:
    """Compute one non-Abelian curvature-scaling summary row."""

    U = su2_midpoint_link_matrix(n_side, potential)
    h = 1.0 / (n_side - 1)
    density = su2_plaquette_log_density(U, n_side)
    her_density = su2_her_reconstructed_log_density(U, n_side)
    target = su2_curvature_values_at_centers(n_side, curvature)
    vector_err = density - target
    norm_err = np.linalg.norm(density, axis=1) - np.linalg.norm(target, axis=1)
    plaquettes = square_plaquettes(n_side, n_side)
    her_errors = np.array([her_plaquette_reconstruction_error(U, p, "su2") for p in plaquettes])
    her_density_err = her_density - density
    log_norms = np.linalg.norm(su2_log_vectors(plaquette_holonomies(U, plaquettes, "su2")), axis=1)
    return SU2CurvatureScalingRow(
        n_side=n_side,
        h=float(h),
        n_plaquettes=len(plaquettes),
        rms_vector_curvature_error=float(np.sqrt(np.mean(np.sum(vector_err**2, axis=1)))),
        max_vector_curvature_error=float(np.max(np.linalg.norm(vector_err, axis=1))),
        rms_norm_curvature_error=float(np.sqrt(np.mean(norm_err**2))),
        max_norm_curvature_error=float(np.max(np.abs(norm_err))),
        max_her_reconstruction_error=float(max(np.max(her_errors), np.max(np.linalg.norm(her_density_err, axis=1)) * h * h)),
        max_plaquette_log_norm=float(np.max(log_norms)) if log_norms.size else 0.0,
    )


def su2_curvature_scaling_table(
    n_sides: Iterable[int] = (8, 12, 16, 24, 32),
    potential: SU2Potential = default_su2_vector_potential,
    curvature: Vector3Field = default_su2_curvature_vector,
) -> list[SU2CurvatureScalingRow]:
    """Compute SU(2) curvature-scaling rows for several grid resolutions."""

    return [su2_curvature_scaling_row(int(n), potential, curvature) for n in n_sides]


def estimate_loglog_slope(h_values: Sequence[float], errors: Sequence[float]) -> float:
    """Estimate slope of log(error) versus log(h)."""

    h = np.asarray(h_values, dtype=float)
    e = np.asarray(errors, dtype=float)
    mask = (h > 0) & (e > 0)
    if int(np.sum(mask)) < 2:
        raise ValueError("Need at least two positive h/error pairs to estimate a slope.")
    slope, _ = np.polyfit(np.log(h[mask]), np.log(e[mask]), deg=1)
    return float(slope)


def max_su2_plaquette_covariance_error_for_smooth_field(
    n_side: int,
    seed: int = 123,
    potential: SU2Potential = default_su2_vector_potential,
) -> float:
    """Check W'_abcd = H_a^-1 W_abcd H_a for a smooth synthetic SU(2) field."""

    U = su2_midpoint_link_matrix(n_side, potential)
    H = random_su2_stack(n_side * n_side, seed=seed)
    return su2_plaquette_covariance_error(U, H, n_side, n_side)


def max_su2_plaquette_trace_invariance_error_for_smooth_field(
    n_side: int,
    seed: int = 123,
    potential: SU2Potential = default_su2_vector_potential,
) -> float:
    """Check trace invariance of smooth SU(2) plaquette holonomies."""

    U = su2_midpoint_link_matrix(n_side, potential)
    H = random_su2_stack(n_side * n_side, seed=seed)
    return su2_plaquette_trace_invariance_error(U, H, n_side, n_side)


def max_su2_her_reconstruction_error_for_smooth_field(
    n_side: int,
    potential: SU2Potential = default_su2_vector_potential,
) -> float:
    """Maximum Her plaquette reconstruction error over a smooth SU(2) lattice."""

    U = su2_midpoint_link_matrix(n_side, potential)
    plaquettes = square_plaquettes(n_side, n_side)
    return float(np.max([her_plaquette_reconstruction_error(U, p, "su2") for p in plaquettes]))
