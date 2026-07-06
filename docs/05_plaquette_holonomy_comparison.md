# Experiment 04: Plaquette holonomy comparison

## Purpose

This experiment compares Her triangle residuals against standard lattice-gauge plaquette holonomy.

The goal is **not** to claim that the Her residual is curvature. The goal is to determine the exact bookkeeping relation between:

1. the standard square Wilson-loop / plaquette holonomy,
2. triangular holonomies obtained by splitting the square along a diagonal, and
3. Her triangle residuals.

## Setup

For an oriented square plaquette

```text
a -> b -> c -> d -> a
```

standard lattice gauge theory uses the loop

```text
W_abcd = U_ab U_bc U_cd U_da.
```

If the diagonal `a -> c` is available, the square decomposes into two based triangular loops:

```text
H_abc = U_ab U_bc U_ca
H_acd = U_ac U_cd U_da
W_abcd = H_abc H_acd.
```

Her's triangle residual is

```text
R_ijk = U_ij^-1 U_ik U_jk^-1.
```

This is not literally the same as `H_ijk`. In the non-Abelian case the relation is:

```text
R_ijk = U_ij^-1 H_ijk^-1 U_ij.
```

So `R_ijk` is conjugate to the inverse triangular holonomy. Equivalently,

```text
H_ijk = U_ij R_ijk^-1 U_ij^-1.
```

## Concrete result

The experiment verifies, for both U(1) and SU(2), that:

```text
W_abcd = H_abc H_acd
```

and that the same plaquette holonomy can be reconstructed from Her residuals by:

```text
W_abcd = (U_ab R_abc^-1 U_ab^-1)(U_ac R_acd^-1 U_ac^-1).
```

## Interpretation

This makes the conservative statement sharper:

> Her triangle residuals are not identical to plaquette curvature. They are basepoint-shifted, inverted triangular consistency residuals whose conjugacy class matches the triangular holonomy class. When a square plaquette is split into two triangles, Her residuals can reconstruct the standard plaquette holonomy exactly, provided the diagonal link is available.

This supports the phrase **plaquette-compatible relational consistency diagnostic**.

It does not yet support the phrase **new curvature law**.

## Boundary

The result depends on the availability of the diagonal link `U_ac`. In a nearest-neighbor-only lattice, that diagonal is not automatically present. It must be supplied as an additional link, a composed path, a chart transition, or an experimentally meaningful relational measurement.

That distinction matters. Without the diagonal, the Her triangle residual and the standard square plaquette live on different combinatorial data.
