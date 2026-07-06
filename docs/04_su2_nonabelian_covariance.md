# SU(2) Non-Abelian Covariance Layer

This layer tests whether the U(1) re-anchoring result survives in a genuinely non-Abelian setting.

## Convention

Directed links are written

\[
U_{ij}: i \to j.
\]

A local SU(2) gauge transformation acts by

\[
U_{ij}\mapsto U'_{ij}=H_i^{-1}U_{ij}H_j.
\]

For pure gauge site frames \(G_i\), define

\[
U_{ij}=G_i^{-1}G_j.
\]

Then Her re-anchoring predicts

\[
U_{ij}^{-1}U_{ik}=U_{jk}.
\]

The triangle residual is

\[
R_{ijk}=U_{ij}^{-1}U_{ik}U_{jk}^{-1}.
\]

Pure gauge consistency means

\[
R_{ijk}=I.
\]

## The important non-Abelian distinction

In U(1), residuals are literally invariant under gauge transformation because all phases commute.

In SU(2), residual matrices are not generally literally invariant. Instead,

\[
R'_{ijk}=H_j^{-1}R_{ijk}H_j.
\]

So the matrix residual is **gauge-covariant by conjugation at the re-anchored endpoint** \(j\).

Gauge-invariant diagnostics include:

- \(\operatorname{tr}(R_{ijk})\),
- \(\lVert R_{ijk}-I\rVert_F\),
- eigenvalues of \(R_{ijk}\),
- the SU(2) rotation angle inferred from \(\operatorname{tr}(R_{ijk})\).

## What this supports

This supports the concrete statement:

> Her re-anchoring defines a non-Abelian triangle consistency residual. For SU(2) pure gauges, the residual vanishes to numerical precision. Under local SU(2) gauge transformations, the residual transforms covariantly by endpoint conjugation rather than remaining literally invariant.

## What this does not support yet

This still does not prove that Her Algorithm is a new physical theory. It shows that the identity can be implemented coherently inside a non-Abelian lattice-gauge-style setting.

Before making curvature claims, the next required comparison is with standard plaquette/Wilson-loop holonomy.
