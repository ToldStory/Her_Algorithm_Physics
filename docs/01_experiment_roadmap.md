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
