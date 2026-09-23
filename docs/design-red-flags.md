# Design red flags

Problems found while documenting the code. Each entry: location, problem.
Nothing here has been fixed yet.

- `src/model/mesh/api/include/model.h`, `getModelTheta*` / `getModelPhi*`: return
  `ScalarType` (int) although `ModelUnstruct` stores angles in radians and
  `sem_solver_impl.h` reads them as float. Angles are truncated; TTI is wrong for
  any angle below 1 rad on unstructured meshes.
- `src/model/mesh/api/include/model.h`, `ModelApi(const ModelDataBase<ScalarType,
  FloatType>&)`: template arguments swapped with respect to
  `ModelDataBase<FloatType, ScalarType>`; the constructor is unused.
- `ModelApi` and `FaceConnectivityApi` both expose `getNumberOfFaces`,
  `getGlobalFace`, `getGlobalNodeFromFace` and `isBoundaryFace`, with different
  semantics for `isBoundaryFace`: duplicated interface.
- `src/discretization/fe/`: unfinished migration. `QkHexahedronBase` is never
  derived from, `FeDiscretizationTag` is never inherited, `AssertFeDiscretization`
  and `DiscretizationTraits` are never used; dispatch still goes through the
  deprecated `IntegralTypeSelector`.
- Mass matrix assembly (`computeGlobalMassMatrix` in `sem_solver_impl.h`,
  `initGeometricMassMatrix` in both differentiators): the local index is decoded
  as x = i % n, z = (i/n) % n, y = i/n^2, then passed to globalNodeIndex(e, x, y, z),
  which swaps y and z relative to computeMassTerm's linearIndex(qa, qb, qc).
  Invisible on axis-aligned boxes, suspected wrong on distorted hexahedra.
