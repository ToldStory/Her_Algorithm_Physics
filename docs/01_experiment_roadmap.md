# Experiment Roadmap

1. Reproduce group law.
2. Build \(U(1)\) lattice gauge fields.
3. Test gauge covariance.
4. Test SU(2) non-Abelian covariance by conjugation.
5. Compare residuals to plaquette holonomy.
6. Generate flat versus curved synthetic fields.
7. Localize corrupted links.
8. Test clustered faults.
9. Study non-invertible and near-singular boundaries.


## 04. Plaquette holonomy comparison

Branch:

    investigate/plaquette-holonomy-comparison

Goal:

Compare Her triangle residuals against standard square plaquette holonomy / Wilson-loop diagnostics.

Tests:

- pure-gauge plaquette flatness in U(1) and SU(2),
- exact square factorization into two triangular holonomies,
- exact reconstruction of square holonomy from Her residuals,
- U(1) gauge invariance,
- SU(2) gauge covariance by basepoint conjugation,
- trace / angle invariance as conjugacy diagnostics.

Critical boundary:

This does not prove that the Her residual is curvature. It proves that Her residuals are exactly compatible with plaquette holonomy after diagonal splitting, basepoint transport, and inversion bookkeeping.


## 05. Curvature scaling limit

Goal:

Test whether Her-reconstructed plaquette holonomy obeys the standard small-area U(1) curvature scaling law when links are sampled from a smooth synthetic vector potential.

Tests:

- constant curvature field \(A=(0,Bx)\) produces plaquette angle density \(B\),
- nonconstant smooth field converges under lattice refinement,
- estimated log-log error slope is positive and near second order for midpoint links,
- Her-reconstructed plaquette holonomy matches ordinary plaquette holonomy,
- U(1) plaquette values remain gauge-invariant under local phase changes.

Critical boundary:

This does not prove Her residuals are literally curvature. It proves that, after reconstructing plaquette holonomy from Her residuals, the resulting plaquette diagnostic has the expected smooth U(1) curvature scaling behavior.


## 06. SU(2) curvature scaling limit

Goal:

Test whether Her-reconstructed plaquette holonomy obeys the standard small-area SU(2) curvature scaling law when links are sampled from a smooth synthetic non-Abelian gauge potential.

Tests:

- smooth SU(2) connection represented by Pauli-vector fields,
- midpoint-sampled SU(2) links,
- principal logarithm vectors `log(W_square)/h^2`,
- comparison against analytic non-Abelian curvature including the commutator term,
- estimated convergence under lattice refinement,
- exact Her reconstruction of plaquette holonomy,
- SU(2) plaquette covariance and trace invariance under local gauge transformations.

Critical boundary:

This still does not prove that Her residuals are literally curvature in all settings. It proves that, after reconstructing plaquette holonomy from Her residuals, the resulting non-Abelian plaquette diagnostic has the expected smooth SU(2) small-area curvature behavior in a controlled synthetic test.
