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
python experiments/02_u1_lattice_gauge.py
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


## Second experiment: SU(2) non-Abelian covariance

After locking the U(1) baseline, create and work on:

```bash
git checkout -b investigate/su2-nonabelian-covariance
```

Run:

```bash
python experiments/03_su2_nonabelian_covariance.py
```

This tests:

- pure-gauge SU(2) links \(U_{ij}=G_i^{-1}G_j\),
- non-Abelian re-anchoring \(U_{ij}^{-1}U_{ik}=U_{jk}\),
- residuals \(R_{ijk}=U_{ij}^{-1}U_{ik}U_{jk}^{-1}\),
- covariance \(R'_{ijk}=H_j^{-1}R_{ijk}H_j\),
- failure of literal matrix invariance for generic non-Abelian residuals,
- trace / rotation-angle invariance under conjugation,
- sparse perturbed-link localization by residual participation.

Outputs:

```text
outputs/tables/su2_nonabelian_covariance_summary.csv
outputs/figures/su2_perturbed_residual_histogram.png
```


## Third experiment: plaquette / holonomy comparison

After locking the SU(2) covariance layer, create and work on:

    git checkout -b investigate/plaquette-holonomy-comparison

Run:

    python experiments/04_plaquette_holonomy_comparison.py

This tests:
  * standard square plaquette holonomy W_abcd = U_ab U_bc U_cd U_da,
  * exact diagonal factorization W_abcd = H_abc H_acd,
  * relation between Her residuals and triangular holonomy,
  * reconstruction of plaquette holonomy from Her residuals,
  * U(1) plaquette gauge invariance,
  * SU(2) plaquette gauge covariance and trace invariance.

Outputs:
    outputs/tables/plaquette_holonomy_comparison_summary.csv
    outputs/figures/u1_su2_plaquette_angles.png


## Fourth experiment: curvature scaling limit

Work directly on `main` after confirming the current baseline passes tests.

Run:

```bash
python experiments/05_curvature_scaling_limit.py
```

This tests:

- smooth synthetic U(1) midpoint-sampled links,
- plaquette angle scaling against analytic curvature \(F_{xy}\),
- convergence of \(rg(W_\square)/h^2\) toward \(F_{xy}\) under lattice refinement,
- exact Her reconstruction of the same plaquette holonomy,
- gauge invariance of the U(1) plaquette diagnostic.

Outputs:

```text
outputs/tables/curvature_scaling_limit_summary.csv
outputs/figures/curvature_scaling_limit.png
```

Safe claim:

> Her residuals reconstruct plaquette holonomy, and in this synthetic U(1) refinement test the reconstructed plaquette holonomy follows the expected curvature-times-area scaling.

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
