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
- convergence of \(\arg(W_\square)/h^2\) toward \(F_{xy}\) under lattice refinement,
- exact Her reconstruction of the same plaquette holonomy,
- gauge invariance of the U(1) plaquette diagnostic.

Outputs:

```text
outputs/tables/curvature_scaling_limit_summary.csv
outputs/figures/curvature_scaling_limit.png
```

Safe claim:

> Her residuals reconstruct plaquette holonomy, and in this synthetic U(1) refinement test the reconstructed plaquette holonomy follows the expected curvature-times-area scaling.

## Fifth experiment: SU(2) curvature scaling limit

Work directly on `main` after confirming the current baseline passes tests.

Run:

```bash
python experiments/06_su2_curvature_scaling_limit.py
```

This tests:

- smooth synthetic SU(2) midpoint-sampled links,
- principal logarithm vectors `log(W_square)/h^2` against analytic non-Abelian curvature,
- inclusion of the SU(2) commutator term in the analytic target,
- convergence of the plaquette diagnostic under lattice refinement, including the non-Abelian basepoint subtlety,
- exact Her reconstruction of the same plaquette holonomy,
- SU(2) plaquette covariance and trace invariance under local gauge transformations.

Outputs:

```text
outputs/tables/su2_curvature_scaling_limit_summary.csv
outputs/figures/su2_curvature_scaling_limit.png
```

Safe claim:

> Her residuals reconstruct non-Abelian plaquette holonomy, and in this synthetic SU(2) refinement test the reconstructed plaquette holonomy follows the expected curvature-times-area scaling. The raw Lie-algebra vector comparison converges, while the conjugacy-invariant norm diagnostic shows near-second-order scaling and includes the commutator contribution to curvature.

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

Established Results

This repository investigates Her Algorithm as a discrete relational re-anchoring structure for group-valued fields. The central algebraic identity is exact: for group elements a,b,c∈G, the relational differences

ϕ(a,b)=a
−1
b,ϕ(a,c)=a
−1
c

satisfy

ϕ(a,b)
−1
ϕ(a,c)=b
−1
c.

Thus, a common origin or anchor a cancels without requiring commutativity. This makes the construction naturally compatible with both Abelian and non-Abelian settings.

The work in this repository does not claim that Her Algorithm is a new physical theory. Instead, it establishes a sequence of concrete mathematical and computational facts showing that the re-anchoring identity behaves correctly inside lattice-gauge-style systems.

1. U(1) Gauge Re-Anchoring

The first experiment implements Her re-anchoring in a finite U(1) lattice-gauge setting. For pure gauge fields, the re-anchoring residual vanishes to numerical precision. When sparse synthetic perturbations are introduced, the residual becomes nonzero and localizes the affected relational inconsistencies.

Because U(1) is Abelian, the residual diagnostics are gauge-invariant under local U(1) transformations.

Established statement:

In finite U(1) link fields, Her re-anchoring defines a triangle consistency residual that vanishes for pure gauges, detects synthetic perturbations, and remains invariant under local gauge transformations.

2. SU(2) Non-Abelian Covariance

The second experiment extends the construction to SU(2), where multiplication is noncommutative. This is a stronger test because residual matrices are not expected to remain literally unchanged under local gauge transformations.

The result confirms the expected non-Abelian behavior: the residual transforms covariantly by conjugation, while conjugacy-invariant diagnostics such as trace, Frobenius distance from identity, and rotation angle remain invariant.

Established statement:

In finite SU(2) link fields, Her re-anchoring defines a non-Abelian triangle consistency residual. The residual vanishes for pure gauges and transforms covariantly under local SU(2) gauge transformations. Its conjugacy-invariant observables remain gauge-invariant.

This distinguishes the SU(2) case from the U(1) case. In U(1), residuals are literally invariant. In SU(2), the residual matrices may change substantially, but they change in the mathematically correct way.

3. Plaquette and Holonomy Compatibility

The third experiment compares Her triangle residuals with standard lattice-gauge plaquette holonomy.

For a square plaquette with vertices a,b,c,d, the standard plaquette holonomy is

W
abcd
	​

=U
ab
	​

U
bc
	​

U
cd
	​

U
da
	​

.

When the square is split into two triangles along a diagonal, the corresponding triangular holonomies reconstruct the square plaquette holonomy:

W
abcd
	​

=H
abc
	​

H
acd
	​

.

Her residuals are not identical to plaquette curvature by themselves. Rather, they are basepoint-shifted and inverted triangle consistency residuals. With correct inversion and basepoint bookkeeping, they reconstruct the triangular holonomies, and therefore reconstruct the plaquette holonomy.

Established statement:

Her triangle residuals are plaquette-compatible. When a square plaquette is decomposed into two triangles, Her residuals reconstruct the standard plaquette holonomy up to numerical precision, after inversion and basepoint bookkeeping.

This result is important because it clarifies the relationship between Her re-anchoring and lattice gauge theory. The residual is not merely an arbitrary defect score. It is structurally compatible with Wilson-loop / plaquette holonomy.

4. Curvature Scaling Limit

The fourth experiment tests whether the Her-reconstructed plaquette holonomy behaves like a curvature diagnostic under lattice refinement.

A smooth synthetic U(1) connection is generated on increasingly fine lattices. The plaquette angle divided by the lattice area is compared against the analytic curvature F
xy
	​

. The observed log-log convergence slope is approximately

1.995,

which is essentially second-order convergence.

Established statement:

In a smooth synthetic U(1) connection, the Her-reconstructed plaquette holonomy reproduces the expected curvature-times-area behavior under lattice refinement, with approximately second-order convergence.

This is the strongest result so far. It shows that the construction is not only algebraically compatible with plaquette holonomy, but also behaves correctly in a smooth curvature-scaling test.

Current Claim Boundary

The repository currently supports the following disciplined claim:

Her Algorithm is an exact relational re-anchoring identity for group-valued data. In lattice-gauge-style systems, this identity defines triangle consistency residuals that vanish on pure gauges, detect perturbations, transform correctly under local gauge transformations, reconstruct plaquette holonomy, and reproduce expected curvature-times-area scaling in a smooth U(1) refinement test.

The repository does not yet claim:

Her Algorithm is a new physical law.

Nor does it claim:

Her residuals are literally curvature in all settings.

The safer and currently supported statement is:

Her re-anchoring is a gauge-compatible relational diagnostic. Through plaquette reconstruction, it can reproduce standard lattice-gauge holonomy behavior and, in the tested U(1) case, the expected curvature scaling limit.

Verified Test State

The current codebase contains tests for:

U(1) gauge re-anchoring
SU(2) non-Abelian covariance
Plaquette / holonomy reconstruction
Curvature scaling limit

The current validated test state is:

25 passed

This means the project has moved beyond informal analogy. It now contains a reproducible computational sequence supporting precise claims about relational re-anchoring, gauge covariance, plaquette compatibility, and curvature scaling behavior.
