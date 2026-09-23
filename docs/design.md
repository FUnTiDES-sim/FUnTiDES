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

## 1D Lagrange bases (`LagrangeBasis*`)

`LagrangeBasis1`, `LagrangeBasis2` and `LagrangeBasisNGL` (N = 3 to 9) are stateless classes with
the same static interface, used as the `GL_BASIS` template argument of the hexahedral
discretizations. For a basis of order `r`, with `n = r + 1` nodes:

- The parent interval is `[-1, 1]`. Its nodes are the `n` Gauss-Lobatto-Legendre points, in
  increasing order: node 0 is -1, node `r` is +1 (for `r` = 1 and 2 they are equispaced).
  `parentSupportCoord(q)` returns node `q`, `weight(q)` its quadrature weight (the weights sum
  to 2).
- `value(q, xi)` is the Lagrange polynomial of node `q` (1 at node `q`, 0 at the other nodes);
  `gradient(q, xi)` is its derivative with respect to `xi`. `valueK()` / `gradientK()` are the
  same for a fixed `q = K`. Indices must lie in `[0, r]`; no bound is checked.
- `gradientAt(q, p)` is `gradient(q, parentSupportCoord(p))`, precomputed. Callers only query
  `p <= (n - 1) / 2` and obtain the other half from the symmetry
  `gradientAt(q, p) = -gradientAt(r - q, r - p)`. For `r <= 5` the table only holds that half
  and returns meaningless values beyond it; for `r >= 6` it is complete.
- `TensorProduct2D` and `TensorProduct3D` give the tensor-product nodes: `linearIndex(i, j)` is
  `i + n*j` and `linearIndex(i, j, k)` is `i + n*j + n*n*k` (hexahedron local numbering above),
  `multiIndex()` is its inverse, and `value(coords, N)` fills `N[linearIndex]` with the product
  of the 1D basis values at the parent coordinates `coords`.

The coefficients of orders 6 to 9 were produced by a generator script (`computGL.py`) that is
not in the repository.

## Symmetric 3x3 matrices (Voigt storage)

The discretization kernels store a symmetric 3x3 matrix `A` as 6 values in the order
`(A00, A11, A22, A12, A02, A01)`, and a symmetric 2x2 matrix as `(A00, A11, A01)`.
