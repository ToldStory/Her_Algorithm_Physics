"""Plaquette / Wilson-loop comparisons for Her re-anchoring residuals.

This module is intentionally conservative.  It does not define Her residuals
as curvature.  It compares Her triangle residuals against the standard
lattice-gauge plaquette holonomy on a square that has been split into two
triangles by a diagonal.

Conventions
-----------
A directed link U[i, j] transports from site i to site j.

For a square plaquette (a,b,c,d), oriented counterclockwise as

a -> b -> c -> d -> a,

the standard plaquette holonomy is

    W_abcd = U_ab U_bc U_cd U_da.

If the diagonal U_ac is available, the same loop factors exactly into two
based triangular holonomies:

    W_abcd = H_abc H_acd,

where

    H_xyz = U_xy U_yz U_zx.

Her's triangle residual is

    R_ijk = U_ij^-1 U_ik U_jk^-1.

It is not literally the same matrix as H_ijk.  In the non-Abelian case,
R_ijk is conjugate to the inverse triangular holonomy:

    R_ijk = U_ij^-1 H_ijk^-1 U_ij.

Therefore Her residuals preserve the conjugacy class of triangular holonomy,
and, when combined with the relevant link, can reconstruct the based triangular
holonomy exactly:

    H_ijk = U_ij R_ijk^-1 U_ij^-1.

The tests and experiment use this identity to determine whether Her residuals
are a curvature proxy, a consistency diagnostic, or a separate relational
object. The first result is deliberately narrow: Her residuals are exactly
compatible with plaquette holonomy after basepoint/inversion bookkeeping.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np

from .u1_gauge import (
    angle_errors as u1_angle_errors,
    gauge_transform as u1_gauge_transform,
    inverse_u1,
    triangle_residual as u1_triangle_residual,
)
from .su2_gauge import (
    I2,
    dagger,
    gauge_transform as su2_gauge_transform,
    inverse_su2,
    su2_rotation_angles,
    triangle_residual as su2_triangle_residual,
)

Array = np.ndarray
GroupKind = Literal["u1", "su2"]


@dataclass(frozen=True)
class PlaquetteSummary:
    """Summary statistics for a set of oriented plaquettes."""

    group: str
    nx: int
    ny: int
    n_plaquettes: int
    max_plaquette_angle: float
    mean_plaquette_angle: float
    max_her_reconstruction_error: float
    max_triangle_factorization_error: float


def grid_site(x: int, y: int, nx: int) -> int:
    """Row-major site index for a rectangular grid."""

    return y * nx + x


def square_plaquettes(nx: int, ny: int) -> list[tuple[int, int, int, int]]:
    """Return oriented square plaquettes (a,b,c,d) on an nx-by-ny grid.

    The orientation is

        a=(x,y) -> b=(x+1,y) -> c=(x+1,y+1) -> d=(x,y+1) -> a.
    """

    if nx < 2 or ny < 2:
        raise ValueError("A plaquette grid needs nx >= 2 and ny >= 2.")
    plaquettes: list[tuple[int, int, int, int]] = []
    for y in range(ny - 1):
        for x in range(nx - 1):
            a = grid_site(x, y, nx)
            b = grid_site(x + 1, y, nx)
            c = grid_site(x + 1, y + 1, nx)
            d = grid_site(x, y + 1, nx)
            plaquettes.append((a, b, c, d))
    return plaquettes


def plaquette_edges(nx: int, ny: int) -> list[tuple[int, int]]:
    """Return undirected nearest-neighbor grid edges in canonical orientation."""

    edges: list[tuple[int, int]] = []
    for y in range(ny):
        for x in range(nx):
            p = grid_site(x, y, nx)
            if x + 1 < nx:
                edges.append((p, grid_site(x + 1, y, nx)))
            if y + 1 < ny:
                edges.append((p, grid_site(x, y + 1, nx)))
    return edges


def _inv(Uij: Array | complex, group: GroupKind) -> Array | complex:
    if group == "u1":
        return inverse_u1(Uij)  # type: ignore[arg-type]
    if group == "su2":
        return inverse_su2(Uij)  # type: ignore[arg-type]
    raise ValueError(f"Unknown group kind: {group!r}")


def _mul3(A: Array | complex, B: Array | complex, C: Array | complex) -> Array | complex:
    if np.ndim(A) == 0:
        return A * B * C  # type: ignore[operator]
    return A @ B @ C  # type: ignore[operator]


def _mul4(
    A: Array | complex,
    B: Array | complex,
    C: Array | complex,
    D: Array | complex,
) -> Array | complex:
    if np.ndim(A) == 0:
        return A * B * C * D  # type: ignore[operator]
    return A @ B @ C @ D  # type: ignore[operator]


def _identity_like(group: GroupKind) -> Array | complex:
    return 1.0 + 0.0j if group == "u1" else I2


def _norm(A: Array | complex) -> float:
    return float(abs(A)) if np.ndim(A) == 0 else float(np.linalg.norm(A))


def triangular_holonomy(
    U: Array,
    i: int,
    j: int,
    k: int,
    group: GroupKind,
) -> Array | complex:
    """Based triangular holonomy H_ijk = U_ij U_jk U_ki."""

    return _mul3(U[i, j], U[j, k], _inv(U[i, k], group))


def plaquette_holonomy(
    U: Array,
    plaquette: tuple[int, int, int, int],
    group: GroupKind,
) -> Array | complex:
    """Oriented plaquette holonomy W_abcd = U_ab U_bc U_cd U_da."""

    a, b, c, d = plaquette
    return _mul4(U[a, b], U[b, c], U[c, d], _inv(U[a, d], group))


def plaquette_holonomies(
    U: Array,
    plaquettes: Sequence[tuple[int, int, int, int]],
    group: GroupKind,
) -> Array:
    """Compute holonomies for a list of oriented plaquettes."""

    return np.array([plaquette_holonomy(U, p, group) for p in plaquettes])


def split_plaquette_holonomy(
    U: Array,
    plaquette: tuple[int, int, int, int],
    group: GroupKind,
) -> Array | complex:
    """Compute W_abcd as H_abc H_acd using diagonal a-c."""

    a, b, c, d = plaquette
    H1 = triangular_holonomy(U, a, b, c, group)
    H2 = triangular_holonomy(U, a, c, d, group)
    if group == "u1":
        return H1 * H2  # type: ignore[operator]
    return H1 @ H2  # type: ignore[operator]


def triangle_factorization_error(
    U: Array,
    plaquette: tuple[int, int, int, int],
    group: GroupKind,
) -> float:
    """Norm of W_abcd - H_abc H_acd."""

    return _norm(plaquette_holonomy(U, plaquette, group) - split_plaquette_holonomy(U, plaquette, group))


def her_triangle_residual(
    U: Array,
    i: int,
    j: int,
    k: int,
    group: GroupKind,
) -> Array | complex:
    """Dispatch to the U(1) or SU(2) Her triangle residual."""

    if group == "u1":
        return u1_triangle_residual(U, i, j, k)
    if group == "su2":
        return su2_triangle_residual(U, i, j, k)
    raise ValueError(f"Unknown group kind: {group!r}")


def triangular_holonomy_from_her(
    U: Array,
    i: int,
    j: int,
    k: int,
    group: GroupKind,
) -> Array | complex:
    """Reconstruct based triangular holonomy from a Her residual.

    Since R_ijk = U_ij^-1 H_ijk^-1 U_ij, we recover

        H_ijk = U_ij R_ijk^-1 U_ij^-1.
    """

    R = her_triangle_residual(U, i, j, k, group)
    if group == "u1":
        return U[i, j] * inverse_u1(R) * inverse_u1(U[i, j])
    return U[i, j] @ inverse_su2(R) @ inverse_su2(U[i, j])


def her_reconstructed_plaquette_holonomy(
    U: Array,
    plaquette: tuple[int, int, int, int],
    group: GroupKind,
) -> Array | complex:
    """Reconstruct W_abcd using Her residuals on triangles abc and acd."""

    a, b, c, d = plaquette
    H1 = triangular_holonomy_from_her(U, a, b, c, group)
    H2 = triangular_holonomy_from_her(U, a, c, d, group)
    if group == "u1":
        return H1 * H2  # type: ignore[operator]
    return H1 @ H2  # type: ignore[operator]


def her_plaquette_reconstruction_error(
    U: Array,
    plaquette: tuple[int, int, int, int],
    group: GroupKind,
) -> float:
    """Norm of W_abcd minus the Her-reconstructed plaquette holonomy."""

    W = plaquette_holonomy(U, plaquette, group)
    W_her = her_reconstructed_plaquette_holonomy(U, plaquette, group)
    return _norm(W - W_her)


def her_triangle_holonomy_conjugacy_error(
    U: Array,
    i: int,
    j: int,
    k: int,
    group: GroupKind,
) -> float:
    """Check R_ijk = U_ij^-1 H_ijk^-1 U_ij."""

    R = her_triangle_residual(U, i, j, k, group)
    H = triangular_holonomy(U, i, j, k, group)
    predicted = _mul3(_inv(U[i, j], group), _inv(H, group), U[i, j])
    return _norm(R - predicted)


def plaquette_angles(holonomies: Array, group: GroupKind) -> Array:
    """Gauge-invariant angular size of plaquette holonomies."""

    if group == "u1":
        return u1_angle_errors(holonomies)
    if group == "su2":
        return su2_rotation_angles(holonomies)
    raise ValueError(f"Unknown group kind: {group!r}")


def plaquette_summary(
    U: Array,
    nx: int,
    ny: int,
    group: GroupKind,
) -> PlaquetteSummary:
    """Summarize plaquette holonomy and Her reconstruction errors."""

    plaquettes = square_plaquettes(nx, ny)
    W = plaquette_holonomies(U, plaquettes, group)
    angles = plaquette_angles(W, group)
    her_errors = np.array([her_plaquette_reconstruction_error(U, p, group) for p in plaquettes])
    split_errors = np.array([triangle_factorization_error(U, p, group) for p in plaquettes])
    return PlaquetteSummary(
        group=group,
        nx=nx,
        ny=ny,
        n_plaquettes=len(plaquettes),
        max_plaquette_angle=float(np.max(angles)) if angles.size else 0.0,
        mean_plaquette_angle=float(np.mean(angles)) if angles.size else 0.0,
        max_her_reconstruction_error=float(np.max(her_errors)) if her_errors.size else 0.0,
        max_triangle_factorization_error=float(np.max(split_errors)) if split_errors.size else 0.0,
    )


def u1_plaquette_gauge_invariance_error(
    U: Array,
    lambdas: Sequence[float],
    nx: int,
    ny: int,
) -> float:
    """Max literal U(1) plaquette change under local gauge transformation."""

    plaquettes = square_plaquettes(nx, ny)
    W = plaquette_holonomies(U, plaquettes, "u1")
    V = u1_gauge_transform(U, lambdas)
    W2 = plaquette_holonomies(V, plaquettes, "u1")
    return float(np.max(np.abs(W2 - W))) if len(plaquettes) else 0.0


def su2_plaquette_covariance_error(U: Array, H: Array, nx: int, ny: int) -> float:
    """Max error in W'_abcd = H_a^-1 W_abcd H_a for SU(2)."""

    plaquettes = square_plaquettes(nx, ny)
    W = plaquette_holonomies(U, plaquettes, "su2")
    V = su2_gauge_transform(U, H)
    W2 = plaquette_holonomies(V, plaquettes, "su2")
    H_inv = dagger(H)
    errs = []
    for idx, (a, _, _, _) in enumerate(plaquettes):
        predicted = H_inv[a] @ W[idx] @ H[a]
        errs.append(np.linalg.norm(W2[idx] - predicted))
    return float(np.max(errs)) if errs else 0.0


def su2_plaquette_trace_invariance_error(U: Array, H: Array, nx: int, ny: int) -> float:
    """Max |tr(W') - tr(W)| under local SU(2) gauge transformation."""

    plaquettes = square_plaquettes(nx, ny)
    W = plaquette_holonomies(U, plaquettes, "su2")
    V = su2_gauge_transform(U, H)
    W2 = plaquette_holonomies(V, plaquettes, "su2")
    return float(np.max(np.abs(np.trace(W2, axis1=-2, axis2=-1) - np.trace(W, axis1=-2, axis2=-1))))


def plaquette_participation_scores(
    U: Array,
    nx: int,
    ny: int,
    group: GroupKind,
) -> dict[tuple[int, int, int, int], float]:
    """Score each square plaquette by its gauge-invariant holonomy angle."""

    plaquettes = square_plaquettes(nx, ny)
    W = plaquette_holonomies(U, plaquettes, group)
    angles = plaquette_angles(W, group)
    return {p: float(angle) for p, angle in zip(plaquettes, angles)}


def top_scored_plaquettes(
    scores: dict[tuple[int, int, int, int], float],
    m: int,
) -> list[tuple[tuple[int, int, int, int], float]]:
    """Return top-m plaquettes by descending score."""

    return sorted(scores.items(), key=lambda item: item[1], reverse=True)[:m]
