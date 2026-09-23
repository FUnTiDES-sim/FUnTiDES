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
- TTI angle unit: `computeCTensor()` (`elasticity_utils.h`, used by
  `ModelApi::initElasticityTensors()`) converts theta and phi from degrees, while
  `computeCMatrix()` in `sem_solver.h`, fed by `getModelTheta/PhiOnNodes`, documents
  radians. The same model data is read in two units depending on the code path.
- `src/model/mesh/api/include/model.h`, `ModelApi::nodeCoord`: `ModelStruct` adds the
  subdomain origin, `ModelUnstruct` returns coordinates built with a zero offset
  (`CartesianUnstructBuilder::initNodesCoords`), so unstructured subdomains of
  different ranks overlap in space.
- `ModelApi::getMinSpacing`: `ModelUnstruct` only measures element 0, which is wrong
  for non-uniform meshes.
- `ModelStruct::getMaxSpeed` returns a constant 1500 and
  `ModelStruct::initElasticityTensors` builds the TTI tensor from hardcoded
  vp/vs/rho, both ignoring the per-node/per-element arrays the model may hold.
- `BoundaryFlag::Sponge` and `BoundaryFlag::Ghost` are never assigned (only exposed
  to Python): planned architecture not implemented.
- `ModelBuilderBase` and `ModelApi` are polymorphic bases with a non-virtual
  destructor; deleting a derived object through a base pointer is undefined.
- `ModelBuilderBase::MAX_ORDER` and `MAX_GLL_ORDER` (`gllpoints.h`) are two
  independent constants that must stay equal: `CartesianUnstructBuilder` sizes
  buffers with the first and validates the order against the second.
- `src/solver/fe/api/include/solver.h`, `Solver`: `initFEarrays`, `allocateFEarrays`,
  `initSpongeValues`, `resetGlobalVectors`, `computeGlobalMassMatrix` and
  `computeDampingMatrix` are internal steps of `computeFEInit`/`computeForces` exposed as
  public pure virtuals: the interface leaks the SEM implementation.
- `Solver`: the interface is not uniformly implemented. The DG, DG-SEM and p-adaptive
  solvers throw on `getMassMatrix*`, `getDampingMatrix` and `getForceVector`, ignore
  `sponge_size`/`surface_sponge`/`taper_delta` in `computeFEInit` and silently ignore
  `setSLSAttenuation`; each `outputSolutionValues` overload is a no-op in one solver family.
- `Solver::computeFEInit`, `surface_sponge`: contradictory meaning. The CLI help says
  surface nodes are "non sponge nodes", `SemProxy::surface_sponge_` says "the top surface has
  an absorbing boundary", and `SEMsolver::initSpongeValues` drops the sponge on the x = 0
  side (not z).
- `Solver::setZBoundary` and `Solver::setElementTags`: hooks for two coupled solvers on the
  base interface, whose tag values are defined per implementation (`kElementTypeSEM`,
  `kElementTypePMin`...), so a caller cannot use them through `Solver` alone.
- `Solver::getInterfaceCouplingCoeff`: the default returns a non-const reference to a
  function-local static shared by all solvers; a caller can write into it.
- `src/solver/fe/api/include/rhs.h`, `Rhs::getWeights()`: no component index, while elastic
  sources carry one weight set per component (`RhsElastic::getWeights(int)`); the base
  interface cannot express per-component weights.
- `physics_traits.h` and `physics_traits_{acoustic,elastic}.h` exist with the same file
  names in `src/solver/fe` (namespace `solver::fe`) and `src/gradient` (namespace
  `gradient`); which one `#include "physics_traits.h"` picks depends on include path order.
- `src/gradient/api/include/wavefield_view.h` and `gradient.h`: the `WavefieldView` and
  `Gradient` bases are never used polymorphically. The differentiator data classes hold the
  concrete types, the Python bindings only expose `print`, and `getNumFields`,
  `getFieldName`, `getNumGradients`, `getGradientName` and `gradient::PhysicsTraits::kName`
  are never called.
- `src/gradient/api/include/differentiator.h`, `Differentiator::compute`, `dt`: used by the
  acoustic differentiator to build the adjoint second time derivative from three snapshots,
  ignored by the elastic one, whose adjoint view carries a precomputed second derivative.
  The two physics expect different adjoint inputs behind the same interface.
- `Differentiator::initGeometricMassMatrix`: documented as required before `compute()`, but
  the acoustic node-based `compute()` builds it lazily and the elastic `compute()` never uses
  it: acoustic node gradients are divided by the nodal volumes, elastic node gradients are
  not. The lazy build also mutates the object from a `const` method through `const_cast`.
