# Design notes

Cross-module conventions. Code refers to the sections below with `@see docs/design.md`.

## Hexahedron local numbering

Applies to every hexahedral element of order `order`, with `n = order + 1` nodes per axis.

- **Element-local DOF index**: `i + j*n + k*n*n`, where `i`, `j`, `k` in `[0, order]` are the
  local indices along x, y and z (x fastest). This is the second index of
  `ModelUnstruct::global_node_index_` and the output of `model::faceLocalToElemLocal()`.
- **CubicFace**: `value / 2` is the normal axis (0 = x, 1 = y, 2 = z), `value % 2` the side
  (0 = minus, 1 = plus), and `value ^ 1` the opposite face.
- **2D face DOF index**: `u + v*n`, where `(u, v)` are the two tangential local indices in axis
  order: `(j, k)` on x faces, `(i, k)` on y faces, `(i, j)` on z faces. The same ordering is used
  by `model::faceLocalToElemLocal[AtDepth]()`, `FaceConnectivityStruct::getGlobalNodeFromFace()`
  and `FaceConnectivityUnstruct` (which stores face DOFs in the owner element's ordering).
- **Structured Cartesian grid** (`ModelStruct`, `FaceConnectivityStruct`): element
  `e = ei + ej*ex + ek*ex*ey`; global node `ix + iy*nx + iz*nx*ny` with `nx = order*ex + 1`,
  `ny = order*ey + 1`.

## Device calls on mesh objects

`ModelApi` and `FaceConnectivityApi` are abstract classes, but virtual dispatch is unusable on
GPU. Kokkos kernels call `PROXY_HOST_DEVICE` methods on the concrete type (solvers and
`FaceConnectivityUnstruct::build()` are templated on the mesh type), never through a base
pointer. Methods without `PROXY_HOST_DEVICE` (e.g. `initElasticityTensors()`,
`buildFaceConnectivity()`, `setQualityFactors()`, `getMaxSpeed()`) are host only.

The same rule applies to `solver::fe::Wavefield` and `solver::fe::Rhs`: device code reaches
them through the concrete types selected by `solver::fe::PhysicsTraits<PHYSICS>`, never through
the base class.

Likewise, `gradient::WavefieldView` and `gradient::Gradient` declare `PROXY_HOST_DEVICE`
virtual getters, but device code may only call them on the concrete types selected by
`gradient::PhysicsTraits<PHYSICS>`.

## Time levels and the split time step

Solvers use an explicit second-order scheme on three time levels. A wavefield holds, per
component, a *current* buffer (u^n), a *previous* buffer (u^(n-1)) and, in backward (adjoint)
mode only, a *previous-previous* buffer. One time step is:

1. `Solver::computeForces()` zeroes and fills the force vectors (`Solver::getForceVector()`),
   one value per global node and component.
2. In a distributed run, the driver sums the force vectors at partition boundaries.
3. `Solver::updateSolutionForward()` writes u^(n+1) into the *previous* buffer;
   `Solver::updateSolutionBackward()` writes it into the *previous-previous* buffer.
4. `Wavefield::swap()` rotates the buffers so that *current* holds u^(n+1) and *previous* u^n.

`Solver::computeOneStep()` performs steps 1 and 3 for non-distributed runs; step 4 is always
the caller's.

## Source terms

`Rhs` stores the point sources of one run. For source `s`: `getElement()[s]` is the element that
contains it, `getTerm(c)(s, t)` is the amplitude of component `c` at time sample `t`, and
`getWeights()(s, l)` the weight of element-local DOF `l` (hexahedron local numbering above).
