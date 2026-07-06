# Her_Algorithm_Physics

This repository investigates **Her Algorithm** as a discrete relational re-anchoring procedure with possible applications to gauge theory, lattice fields, atlas consistency, and physics-inspired diagnostics.

Core law:

\[
(A^{-1}B)^{-1}(A^{-1}C)=B^{-1}C.
\]

Plain language:

\[
A\to B,\quad A\to C \quad\Rightarrow\quad B\to C.
\]

## Critical posture

The theorem above is exact in groups and groupoids. The physics connection is **not** assumed. This repo tests whether the same re-anchoring law is useful as a finite diagnostic for gauge-like link fields, transition consistency, holonomy, curvature, and corrupted-link localization.

The safe public claim is:

> Her Algorithm is a discrete re-anchoring procedure for relational fields, exact where inverse and composition exist.

## First branch

```bash
git checkout -b investigate/u1-gauge-reanchoring
```

## First experiment: U(1) gauge re-anchoring

Install requirements:

```bash
pip install -r requirements.txt
```

Run tests:

```bash
pytest
```

Run the first experiment:

```bash
PYTHONPATH=src python experiments/02_u1_lattice_gauge.py
```

This tests:

- pure-gauge links \(U_{ij}=e^{i(\alpha_j-\alpha_i)}\),
- triangle re-anchoring \(U_{ij}^{-1}U_{ik}=U_{jk}\),
- residuals \(R_{ijk}=U_{ij}^{-1}U_{ik}U_{jk}^{-1}\),
- invariance of U(1) residuals under local gauge transformation,
- sparse perturbed-link residuals,
- preliminary corrupted-link localization by residual participation.

Outputs:

```text
outputs/tables/u1_gauge_reanchoring_summary.csv
outputs/figures/u1_perturbed_residual_histogram.png
```

## What this repo is

A careful mathematical and computational investigation of:

- group-level re-anchoring,
- groupoid re-anchoring,
- gauge-like link variables,
- chart transition consistency,
- cocycle residuals,
- curvature and holonomy diagnostics,
- corrupted-link localization.

## What this repo is not

This repo does not claim Her Algorithm is already a complete physical theory. It does not replace gauge theory. It is not yet validated on real physics data.
