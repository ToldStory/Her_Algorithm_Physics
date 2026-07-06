# Curvature Scaling Limit

This experiment asks a narrower question than “is Her Algorithm curvature?”

It asks:

> If lattice links are sampled from a smooth synthetic U(1) gauge potential, does the Her-reconstructed plaquette holonomy approach the expected small-area curvature law as the lattice spacing shrinks?

## Setup

For a U(1) connection on the unit square, choose a smooth vector potential

\[
A=(A_x,A_y).
\]

The continuum curvature is

\[
F_{xy}=\partial_x A_y-\partial_y A_x.
\]

For a small square plaquette of side length \(h\), ordinary lattice gauge theory predicts

\[
W_\square \approx \exp(iF_{xy}h^2),
\]

so the wrapped plaquette angle divided by \(h^2\) should converge to \(F_{xy}\) at the plaquette center.

## Her reconstruction

The previous plaquette/holonomy layer established that a square plaquette split into two triangles satisfies

\[
W_{abcd}=H_{abc}H_{acd},
\]

and that each triangular holonomy can be reconstructed from the Her residual by

\[
H_{ijk}=U_{ij}R_{ijk}^{-1}U_{ij}^{-1}.
\]

Therefore, in this synthetic U(1) setting, the Her-reconstructed plaquette should have the same small-area scaling as the ordinary plaquette holonomy.

## Result type

A passing result means:

> Her residuals reconstruct plaquette holonomy, and the reconstructed plaquette holonomy obeys the expected U(1) curvature-times-area scaling under lattice refinement.

It does **not** mean:

> Her residuals are literally curvature.

The correct boundary is:

> Her re-anchoring is compatible with curvature diagnostics through plaquette reconstruction in this smooth U(1) scaling test.

## Run

```bash
python -m pytest
python experiments/05_curvature_scaling_limit.py
```

Outputs:

```text
outputs/tables/curvature_scaling_limit_summary.csv
outputs/figures/curvature_scaling_limit.png
```
