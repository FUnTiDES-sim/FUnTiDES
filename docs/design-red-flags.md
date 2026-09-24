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
  derived from, so `FeDiscretizationTag` reaches no back-end and
  `AssertFeDiscretization` would reject both; it and `DiscretizationTraits` are never
  used; dispatch still goes through the deprecated `IntegralTypeSelector`.
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
- `src/discretization/fe/api/fe_discretization_kind.h`, `DiscretizationKind`: declared in
  namespace `solver::fe` inside the discretization API (whose other header uses
  `discretization::fe::api`), and duplicates the runtime selector `utils::enums::implemType`
  (`sem_enums.h`), which only has `kMakutu`: two unrelated enums name the same back-ends.
- `src/io/api/include/io_controller_base.h`, `IOControllerBase`: no solver or driver uses
  it; snapshots still go through the ADIOS2-based `SemIOController`
  (`src/main/fe/include/sem_io_controller.h`). `BackendKind::kAdios2` has no implementation
  (`makeIOController` throws), and `IOConfig::nt`, `IOConfig::nb_receiver` and
  `HostArrayReal` are read by nothing: planned architecture not implemented.
- `IOControllerBase`: the class comment promised row-major storage of multi-dimensional
  arrays independent of `Layout`, but the interface only takes a flat `HostVectorReal` that
  is written verbatim; `local_dims` is metadata never checked against the view size.
- `src/discretization/fe/impl/common/qk_hexahedron_base.h`, `computeMassTerm`,
  `computeStiffnessTerm`, `computeStiffnessTermSumFact`: they take the vertex coordinates as
  `float const (&X)[8][3]` and pass them to `jacobianTransformation` / `computeBMatrix`, which
  expect `real_t const (&)[8][3]`: this cannot compile when `real_t` is `double`. Hidden
  because the class is never instantiated.
- `src/discretization/fe/impl/common/`: the headers are not self-contained. `mathUtilites.h`
  uses `PROXY_HOST_DEVICE` without including `common_macros.h` and mixes it with a second
  host/device macro, `SEMKERNELS_HOST_DEVICE` (`macros.h`); `LagrangeBasis*.h` use `real_t`,
  `PROXY_HOST_DEVICE` and `pow` without including them; `qk_hexahedron_base.h` uses
  `triple_loop` / `for_constexpr`, defined in the makutu back-end header. Conversely
  `Integrals.h` includes both back-end headers, which include `common`: the shared layer
  depends on the back-ends.
- `LagrangeBasis*::gradientAt`: orders 1 to 5 only tabulate the nodes `p <= (n-1)/2` and return
  meaningless values beyond, orders 6 to 9 tabulate every node. The contract depends on the
  order and nothing checks that callers stay in the valid half.
- `src/discretization/fe/impl/common/mathUtilites.h`: `determinant(T const&)`,
  `linearIndex<ORDER>`, `tripleIndex<ORDER>`, `invert3x3` and `computeB` are used nowhere, and
  `src/discretization/fe/impl/makutu/include/tensorops.h`, included by nothing, redefines
  `invert3x3`, `symDeterminant` and `symInvert` with different signatures: dead duplicate code. Its
  `symDeterminant<2>` / `<3>` definitions are partial specializations of function templates,
  which would not compile if the header were included.
- `src/discretization/fe/impl/makutu/include/Qk_Hexahedron_Lagrange_GaussLobatto.h`,
  `computeInterfaceFluxTermAt`: the volume inverse Jacobian is evaluated at the parent point
  `(qa, qb, kQFixed)` for every face, while the face point is `(kQFixed, qa, qb)` on x faces and
  `(qa, kQFixed, qb)` on y faces. Invisible on affine elements (the unit tests only use the unit
  cube), suspected wrong on distorted hexahedra.
- `Qk_Hexahedron_Lagrange_GaussLobatto`: `computeMassTerm`, `computeStiffnessTerm`,
  `computeStiffnessTermSumFact`, `computeStiffNessTermwithJac`, `computeElasticStiffnessSumFact`
  and the vertex overload of `computeElasticStiffnessSumFactTeam` take `float const (&X)[8][3]`,
  and `JacobianType::data` is `float`, but pass them to functions expecting `real_t`: these
  kernels cannot compile when `data_type.h` selects `real_t = double`.
- `Qk_Hexahedron_Lagrange_GaussLobatto`: the mass term uses `|det J|`, while the stiffness kernels
  (`computeBMatrix`, `computeGradPhiGradPhi`, `computeElasticStiffnessSumFact*`) use the signed
  `det J`. For an element whose vertex order gives `det J < 0`, the stiffness changes sign while
  the mass stays positive.
- `Qk_Hexahedron_Lagrange_GaussLobatto::calcGradN`: the quadrature-point overload builds the
  geometry from the 8 vertex nodes of `X[numNodes][3]` only, the parent-coordinate overload from
  all nodes. Same argument, two geometric models.
- `Qk_Hexahedron_Lagrange_GaussLobatto`: dead or legacy members. `getNumQuadraturePoints`,
  `getNumSupportPoints` (non-const) and `getMaxSupportPoints` are virtual in a class with no base,
  a non-virtual destructor and device use; `parentLength` / `parentVolume` are unused and hold
  the first node spacing, not the parent length (2); the header ends with
  `#undef PARENT_GRADIENT_METHOD`, a macro defined nowhere.
- `src/acquisition/include/source_and_receiver_utils.h`, `ComputeRHSWeights`: declares `float invJ[3][3]` and passes it by reference to `invJacobianTransformation`, while the sibling function `ComputeDASWeightsForSample` uses `real_t invJ[3][3]` for the same call; if `real_t` is `double` this reference binding cannot compile, mirroring the known `float`/`real_t` mismatch bug in the Qk_Hexahedron_Lagrange_GaussLobatto kernels.
- `src/acquisition/include/source_and_receiver_utils.h`, `ComputeRHSWeights` vs `ComputeDASWeightsForSample`: two functions computing interpolation weights for a point inside the same element type expose inconsistent interfaces (2D `host_mirror_type` output written to a hardcoded row 0 and overwritten, versus a flat `float*` output that is accumulated), suggesting divergent or duplicated design rather than a single point-evaluation primitive.
- `src/acquisition/include/source_and_receiver_utils.h`: `using namespace std::chrono;` at file scope pollutes every includer's namespace; nothing in the file uses `std::chrono`.
- `src/acquisition/include/source_time_function.h`, `SourceTimeFunction::evaluateRicker`: the default-case message says the order must be 0, 1 or 2 while orders 0 to 4 are implemented, and an unsupported order silently returns 0 after printing to stdout instead of failing.
- `SourceTimeFunction::evaluateRicker`: order 0 uses exp(-2*lam*t^2) while orders 1 to 4 use exp(-lam*t^2), so it is not consistent with the other orders (different Gaussian width); the sign and normalization conventions across orders are also inconsistent.
- `src/acquisition/include/source_time_function.h`: `using namespace std::chrono;` at file scope pollutes every includer and nothing in the file uses `std::chrono`; the header also uses `std::cout`, `std::vector`, `M_PI` and `exp` without including their headers, relying on `data_type.h`.
- `SourceTimeFunction`: the two methods are stateless but neither `static` nor `const`; the cutoff window is computed with double literals (-0.9, 2.9) and the pulse is truncated abruptly, which introduces a discontinuity when tpeak is small relative to 1/f0.
- `src/core/include/common_macros.h`, `FIND_MAX_1D`: expands to an unbraced `if` followed by a second statement, so it is unsafe as the body of an `if`/`else` or loop; it also uses `std::runtime_error` without including `<stdexcept>`, and `FIND_MIN` has no equivalent empty-array check.
- `src/core/include/common_macros.h`, `LOOPHEAD`, `MAINLOOPHEAD`, `KOKKOSNAME`: the macros expand to unbalanced parentheses or a trailing comma and must be paired with `LOOPEND` / `MAINLOOPEND` or placed inside a call; `LaunchMaxThreadsPerBlock` and `LaunchMinBlocksPerSM` are unprefixed macros that pollute the global namespace, and the `LOOPHEAD`/`MAINLOOPHEAD` loop macros are not header-guarded against use outside a class (they capture `this`).
- `src/core/include/data_type_kokkos.h`, `vectorReal`, `arrayReal`, `array3DReal`: despite the "Real" name they are hardcoded to `float` and do not follow `real_t` from `data_type.h`, so they mismatch when `real_t` is `double`.
- `src/core/include/sem_macros.h`, `FENCE`: the expansion already ends with a semicolon, so `FENCE;` yields an empty statement and `if (c) FENCE; else ...` does not compile.
- `src/core/include/sem_macros.h`, `ATOMICADD`: `ADD1` and `ADD2` are not parenthesized in the expansion, so an expression argument may be mis-parsed.
- `src/core/include/sem_macros.h`, `DIMENSION`, `ROW`, `COL`, `ZEROED2D`: unprefixed, generic macro names pollute the global namespace of every includer, and the meaning of `ROW`, `COL` and `ZEROED2D` is not derivable from the header.
- `src/discretization/fe/impl/common/mathUtilites.h`, `symInvert`: the determinant reciprocal is computed as `1.0 / det`, a double expression, so for `T = float` the division is done in double and then narrowed, unlike `invert3x3` which uses `T(1) / det`. Inconsistent precision between the two inverters.
- `src/discretization/fe/impl/common/mathUtilites.h`, `symDeterminant`, `symInvert`, `computeB`: the header defines Voigt-storage helpers without stating the component order in code; the order (0=xx, 1=yy, 2=zz, 3=yz, 4=xz, 5=xy) is only inferable from the formulas.
- `src/discretization/fe/impl/common/qk_hexahedron_base.h`, `interpolationCoord`: it is `constexpr` but not marked `PROXY_HOST_DEVICE`, while it is called from the device-marked `jacobianCoefficient1D`; this relies on implicit constexpr host/device behavior of the compiler.
- `QkHexahedronBase::computeStiffnessTermSumFact`: the declaration names the nodal arrays `p_local` / `f_local` while the definition names them `u_local` / `v_local`; harmless but inconsistent, and the `@param` names follow the declaration.
- `src/discretization/fe/impl/makutu/include/Qk_Hexahedron_Lagrange_GaussLobatto.h`, `interpolationCoord`: `constexpr` but not marked `PROXY_HOST_DEVICE`, while it is called from the device-marked `jacobianCoefficient1D`; relies on implicit constexpr host/device behavior of the compiler.
- `Qk_Hexahedron_Lagrange_GaussLobatto::computeStiffnessTermSumFact`: the declaration names the nodal arrays `p_local` / `f_local` while the definition names them `u_local` / `v_local`; the `@param` names follow the declaration.
- `Qk_Hexahedron_Lagrange_GaussLobatto::computeElasticStiffnessSumFact` and the two `computeElasticStiffnessSumFactTeam` overloads duplicate the same three-pass algorithm, and the single-thread and team versions differ in semantics (accumulate versus overwrite of `f_local`), which is easy to misuse.
- `Qk_Hexahedron_Lagrange_GaussLobatto::computeElasticStiffnessSumFactTeam` (geom overload): the layout of `geom` (10 entries, J^-1 then det J) is a raw-pointer convention with no type or size check, and `F` and `u_local` are unchecked raw pointers.
- `src/discretization/fe/impl/tensorial/include/Qk_Hexahedron_Tensorial.h`, `computeMassTerm`, `computeStiffnessTerm`, `computeStiffnessTermSumFact`, `computeElementMetrics`, `JacobianType::data`: vertex coordinates are `float const (&X)[8][3]` (and `JacobianType` holds `float`) but are passed to `computeBMatrix` / `jacobianTransformation`, which expect `real_t const (&)[8][3]`; this cannot compile when `real_t` is `double`.
- `Qk_Hexahedron_Tensorial_GEMM::computeMassTerm` vs the stiffness kernels: the mass term uses `|det J|` while `computeBMatrix` uses the signed `det J` (via `1/detJ` and `symInvert`), so for a negative Jacobian the stiffness flips sign and the mass does not.
- `Qk_Hexahedron_Tensorial_GEMM`: virtual member functions (`getNumQuadraturePoints`, `getNumSupportPoints`, `getMaxSupportPoints`) in a `final` class with no base and a non-virtual destructor, used in host/device code; the virtual hooks serve no purpose and add a vtable to a device-usable class.
- `Qk_Hexahedron_Tensorial_GEMM::computeDampingTerm`: calls unqualified `sqrt` and `std::abs` without including `<cmath>`, relying on transitive includes.
- `Qk_Hexahedron_Tensorial_GEMM::computeStiffnessOperatorTeamVector`, `computeStiffnessOperatorTeamVectorStreaming`: the input `u` is a non-const `real_t *` although it is only read, and neither function checks the raw-pointer sizes of `u`, `Y`, `W` and `D_flat`.
- `Qk_Hexahedron_Tensorial_GEMM::interpolationCoord`: `constexpr` but not marked `PROXY_HOST_DEVICE` while called from the device-marked `jacobianCoefficient1D`; relies on implicit constexpr host/device behavior.
- `src/gradient/impl/acoustic/include/differentiator_acoustic_impl.h`, `DifferentiatorAcoustic::computeOnElements`, `computeOnNodes`: `invDt2` is declared twice (outside the kernel and again inside the lambda), the outer one being shadowed and unused in `computeOnElements`; the `elementNumber >= getNumberOfElements()` guard inside the `RangePolicy` lambdas is redundant.
- `DifferentiatorAcoustic::initGeometricMassMatrix`, `computeOnElements`, `computeOnNodes`: the kernels hardcode `float` for coordinates and local buffers while the integral callbacks receive `real_t`, so they cannot compile as-is with `real_t = double` (same float/real_t mismatch as the integral back-ends).
- `DifferentiatorAcoustic::computeOnElements`, `computeOnNodes`: gradients are accumulated (`+=` / ATOMICADD) into the output vectors, which are never zeroed here, so the caller must reset them; the accumulate-versus-overwrite contract is not stated by the interface.
- `DifferentiatorAcoustic::initGeometricMassMatrix`: the local buffer is sized `kPointsPerElement` but the loop bound is `mesh.getNumberOfPointsPerElement()`, a runtime value; nothing ties the two together.
- `src/gradient/impl/acoustic/include/gradient_acoustic.h`, `GradientAcoustic::getGradientName`, `getGradient`: an out-of-range index silently returns the kappa name/vector instead of failing, so a wrong index goes undetected.
- `src/gradient/impl/acoustic/include/wavefield_view_backward_acoustic.h`, `WavefieldViewBackwardAcoustic::getField`, `getFieldName`: an out-of-range index silently returns qn (or its name) instead of failing, so a wrong index goes undetected.
- `src/gradient/impl/acoustic/include/wavefield_view_backward_acoustic.h`: the include guard is named `FUNTIDES_GRADIENT_API_INCLUDE_...` although the file lives under `impl/acoustic/include`; the header also uses `vectorReal` and `PROXY_HOST_DEVICE` without including their headers, relying on `wavefield_view.h`.
- `src/gradient/impl/acoustic/include/wavefield_view_forward_acoustic.h`, `WavefieldViewForwardAcoustic::getField`, `getFieldName`: the index is ignored, so any index silently returns pn (or its name) instead of failing.
- `src/gradient/impl/acoustic/include/wavefield_view_forward_acoustic.h`: the include guard is named `FUNTIDES_GRADIENT_API_INCLUDE_...` although the file lives under `impl/acoustic/include`; the header also uses `vectorReal` and `PROXY_HOST_DEVICE` without including their headers, relying on `wavefield_view.h`.
- `WavefieldViewForwardAcoustic`: the constructor is not `explicit`, so a `vectorReal` converts implicitly to a view.
- `src/gradient/impl/common/include/differentiator_factory.h`: the header uses `std::unique_ptr` without including `<memory>`, relying on `differentiator.h`. Its guard and location (`impl/common`) also expose a public factory outside any `api/` directory.
- `src/gradient/impl/common/src/differentiator_factory.cc`, `makeDifferentiatorStruct`, `makeDifferentiatorUnstruct`: the two functions are identical except for `MeshT`, and any physics other than `kAcoustic` silently selects the elastic differentiator (the `// kElastic` else branch has no check).
- `src/gradient/impl/common/src/differentiator_factory.cc`, `orderDispatch`, `makeDifferentiator*`: these templates have external linkage but are not in an anonymous namespace nor `static`, and are not declared in the header, so they are private helpers exposed to the linker. The `MAX_DIFFERENTIATOR_*_ORDER` fallback macros are also defined inside a function body and silently default to 3 when the build definition is missing.
- `src/gradient/impl/common/src/differentiator_factory.cc`, `makeDifferentiatorStruct`, `makeDifferentiatorUnstruct`: `MeshT` hardcodes `float, int` instead of following `real_t`, so the factory cannot follow a `real_t = double` build.
- `src/gradient/impl/common/src/differentiator_factory.cc`, `createDifferentiator`: the parameters `implemType` and `physicType` have the same names as the types they use, which shadows them (`case implemType::kMakutu` works only because of scoping rules for qualified lookup).
- `src/gradient/impl/elastic/include/differentiator_elastic.h`, `DifferentiatorElastic::computeDisplacementGradient`, `computeOnElements`, `computeOnNodes`: coordinates, Jacobians and buffers are hardcoded `float` while the integral back-ends use `real_t`, so the class cannot follow a `real_t = double` build (same float/real_t mismatch as the acoustic differentiator).
- `DifferentiatorElastic::computeOnElements`, `computeOnNodes`: 12 positional `vectorReal` arguments (3 forward, 3 adjoint, 3 second derivative, 3 gradients) with no grouping type; easy to permute silently. The accumulate-versus-overwrite contract for the gradients is not visible from the interface.
- `DifferentiatorElastic::computeOnElements`, `computeOnNodes`: they are public and take `MESH_TYPE` by value, although they are implementation steps of `compute()`.
- `differentiator_elastic.h`, `UNSTRUCT_MESH_TYPE(ORDER)`: ignores its `ORDER` argument, so the same `ModelUnstruct<float, int>` is declared for the three orders; the extern declarations differ only by `ORDER` and `INTEGRAL_TYPE`.
- `src/gradient/impl/elastic/include/differentiator_elastic_impl.h`, `DifferentiatorElastic::computeOnElements`, `computeOnNodes`: the two kernels duplicate about 90 lines (local gather, vertex gather, strain lambda) and differ only in the final scatter; the `dt` parameter is unused in both.
- `DifferentiatorElastic::initGeometricMassMatrix`: `auto mesh = dynamic_cast<MESH_TYPE&>(meshApi)` copies the mesh object by value (`auto` drops the reference); a wrong mesh type throws `std::bad_cast` with no message.
- `DifferentiatorElastic::compute`: `dynamic_cast` to reference on both `data` and `mesh` throws `std::bad_cast` on mismatch, and the six-field backward layout (indices 0-2 adjoint, 3-5 second time derivative) is an implicit convention of `getBackwardField`.
- `src/gradient/impl/elastic/include/gradient_elastic.h`, `GradientElastic::getGradientName`, `getGradient`: an out-of-range index silently returns the rho name/vector instead of failing, so a wrong index goes undetected.
- `src/gradient/impl/elastic/include/gradient_elastic.h`, `GradientElastic`: the header uses `vectorReal` and `PROXY_HOST_DEVICE` without including their headers, relying on `gradient.h`.
- `src/gradient/impl/elastic/include/wavefield_view_backward_elastic.h`, `WavefieldViewBackwardElastic::getField`, `getFieldName`: an out-of-range index silently returns ux_n (or its name) instead of failing, so a wrong index goes undetected.
- `src/gradient/impl/elastic/include/wavefield_view_backward_elastic.h`: the include guard is named `FUNTIDES_GRADIENT_API_INCLUDE_...` although the file lives under `impl/elastic/include`; the header also uses `vectorReal` and `PROXY_HOST_DEVICE` without including their headers, relying on `wavefield_view.h`.
- `WavefieldViewBackwardElastic`: the constructor is not `explicit` (harmless with six arguments), and it takes six positional `vectorReal` arguments that are easy to permute silently.
- `src/gradient/impl/elastic/include/wavefield_view_forward_elastic.h`, `WavefieldViewForwardElastic::getField`, `getFieldName`: an out-of-range index silently returns ux_n (or its name) instead of failing, so a wrong index goes undetected.
- `src/gradient/impl/elastic/include/wavefield_view_forward_elastic.h`: the include guard is named `FUNTIDES_GRADIENT_API_INCLUDE_...` although the file lives under `impl/elastic/include`; the header also uses `vectorReal` and `PROXY_HOST_DEVICE` without including their headers, relying on `wavefield_view.h`.
- `WavefieldViewForwardElastic`: the constructor takes three positional `vectorReal` arguments that are easy to permute silently.
- `src/gradient/pywrap/include/bindings_differentiator.h`, `bind_data_struct`, `bind_gradient_data_acoustic`, `bind_gradient_data_elastic`, `bind_differentiator_base`, `bind_differentiator_factory`: non-template, non-inline functions defined in a header, so including it in more than one translation unit violates the one-definition rule (duplicate symbols at link time). They should be `inline` or moved to a .cc file.
- `bind_differentiator_base`, `get_geometric_mass_matrix`: the lambda returns the Python view by value while the binding uses `reference_internal`; the lifetime tie to the differentiator is only as good as pybind11's handling of this return type, and the registration order (`DataStruct` before derived classes, enums before the factory) is an unchecked implicit requirement.
- `src/core/include/data_type.h`, `real_t`: the alias is declared unconditionally as `float` before the `USE_DOUBLE` selection that redeclares it, so the first declaration is redundant and the header holds three alias declarations for one type.
- `src/core/include/data_type.h`: `using namespace std;` at file scope pollutes every includer, and the header uses unqualified `printf` without including `<cstdio>`.
- `src/core/include/data_type.h`, `allocateVector(int, const char *)`: the `name` parameter is unused (the label is not forwarded to the view), unlike `allocateArray2D` which at least prints it.
- `src/core/include/data_type.h`, `timewatch`, `accumtime`: unprefixed lowercase macros that expand to full statements ending in a semicolon, and `accumtime` accumulates raw clock ticks so the unit depends on the standard library.
- `src/core/include/data_type.h`, `printJMatrix`, `printBMatrix`: debug printing helpers with hardcoded `%f` and an `element < 2` filter live in a core type header (the removed comment itself said they should move elsewhere).
- `src/gradient/pywrap/include/bindings_gradient.h`, `bind_gradient_base`, `bind_gradient_acoustic`, `bind_gradient_elastic`: non-template, non-inline functions defined in a header, so including it in more than one translation unit violates the one-definition rule.
- `bind_gradient_acoustic`, `bind_gradient_elastic`: registering a derived class before `bind_gradient_base` fails at import time; this order is an unchecked implicit requirement.
- `bind_gradient_acoustic`, `bind_gradient_elastic`: the constructors are bound with `python_view_type_t<vectorReal>`, but the gradient classes take `vectorReal`; the Python-side type and the conversion path are not visible from the interface.
- `src/gradient/pywrap/include/bindings_wavefield_view.h`, `bind_wavefield_view_base`, `bind_wavefield_view_forward_acoustic`, `bind_wavefield_view_backward_acoustic`, `bind_wavefield_view_forward_elastic`, `bind_wavefield_view_backward_elastic`: non-template, non-inline functions defined in a header, so including it in more than one translation unit violates the one-definition rule.
- `bind_wavefield_view_*` (derived classes): registering a derived view before `bind_wavefield_view_base` fails at import time; this order is an unchecked implicit requirement.
- `bind_wavefield_view_*` (derived classes): constructors are bound with `python_view_type_t<vectorReal>` while the C++ constructors take `vectorReal`; the conversion path is not visible from the interface. Python users also get only `print`, so the views cannot be used for anything else.
- `src/io/impl/include/posix_io_controller.h`, `PosixIOController::readSnapshot`: takes the destination as `const HostVectorReal&`, although it must write into it (only works because Kokkos views have shallow constness).
- `src/io/impl/include/posix_io_controller.h`: `PosixIOController` is the only working backend but is not used by any driver, and the header sits under `impl/include`.
- `src/io/impl/src/io_controller_factory.cc`, `makeIOController`: the file uses `std::unique_ptr` and `std::make_unique` without including `<memory>`, relying on transitive includes. The public factory is also defined under `impl/src`, outside any `api/` directory.
- `src/io/impl/src/posix_io_controller.cc`, `writeSnapshot`, `readSnapshot`, `makeHeader`, `checkHeader`: the payload size is hardcoded to `sizeof(float)` while `HostVectorReal` follows `real_t`; with `real_t = double` the write truncates the data to half and the read overruns or under-fills the view.
- `PosixIOController::snapshotPath`: the file name is formatted into a fixed `char buf[64]`; a long `prefix` is silently truncated by `snprintf`, so the path loses its index and rank suffix and distinct snapshots can collide.
- `PosixIOController::readSnapshot`: only `nelem` is checked against the view; the stored `dims`, `global_dims`, `offsets` and `snapshot_index` are never compared with the configuration, so a file from a different decomposition is accepted. The header is also written in host byte order with no endianness marker.
- `PosixIOController::flush`: the method is an empty no-op, and `close()` is idempotent only through `closed_`; the interface offers a flush that guarantees nothing.
- `src/main/fe/include/sem_io_controller.h`, `SemIOController::saveReceiver`: `receivers_coords_` is defined with shape {nb_receiver, 3} but each call puts a single `std::array<float, 3>`, so the data written does not match the declared shape unless there is one receiver; the receiver coordinates of the other rows are never provided.
- `SemIOController`: `iter_times_` is defined but never written, and `compressor_op_`, `receiver_op_` and `attachOperator()` are unused placeholders.
- `SemIOController`: the variable names ("AccousticReceiver", "PressureField") and the IO names are hardcoded, misspelled ("Accoustic") and acoustic-specific, so any other field (for example elastic displacement) is stored under the pressure name; element type is hardcoded to `float` while `vectorReal` follows `real_t`.
- `src/main/fe/include/sem_io_controller.h`: the header uses `std::vector`, `std::array` and `std::string` without including `<vector>`, `<array>`, `<string>`, relying on transitive includes. `RECEIVERS_FILE` and `SNAPS_FILE` are unprefixed macros used only to assign file names, and the file names are not configurable.
- `SemIOController`: the class is copyable by default although it owns open ADIOS2 engines closed in the destructor, and the `saveReceiver` receiver argument is a non-const reference although it is only read.
- `src/main/fe/include/sem_proxy.h`, `SEMproxy`: one class mixes acoustic, elastic, acousto-elastic, DG, DG-SEM and p-adaptive state (about 100 members, mostly device arrays and their host mirrors) and orchestration. The interface cannot be summarized in three sentences.
- `SEMproxy::GetPhysic`: returns `int` while its siblings return typed enums (`implemType`, `methodType`, `meshType`). All the `Get*` translators take `std::string` by value.
- `SEMproxy::InitFiniteElem`: the inline body ends with a stray `;` after the closing brace, and it must be called before `Run()` although the constructor already runs setup (`SetupSolver`, `InitMpi`...): two-phase initialization with an unchecked order.
- `SEMproxy::SaveSnapshot`, `snapshot_futures_`: `SaveSnapshot` is `const` but must register a background task, which requires mutating `snapshot_futures_` or a shared state (not `mutable` here). The `const` contract does not match what the method does. The header also includes `<chrono>` without using it directly and receives `using namespace std::chrono` through `source_and_receiver_utils.h` (known).
- `SEMproxy::das_direction_`, `das_vector_`: two 3-vectors describing the same fiber direction (unit and scaled) with no invariant tying them together.
- `SEMproxy::~SEMproxy`: the explicit `io_ctrl_.reset()` is redundant with member destruction; it does not call `WaitSnapshots()`, so pending asynchronous writes are only awaited by the futures' own destructors, which is not guaranteed for futures not created by `std::async`.
- `SEMproxy::dg_sem_iface_z_`, `dg_padaptive_iface_z_`: both default to a hardcoded 1000.f, independent of the domain size.
- `src/main/fe/include/sem_proxy_options.h`, `SemProxyOptions::bind_cli`: several help strings contradict the field defaults. `dt` help says 0.001 s (field 0.006), `timemax` says 1.5 s (field 0.7), `snap-interval` says 10 (field 20), `free-surface` says "Default: true" (field false). The `--implem` help lists only makutu.
- `SemProxyOptions::validate`: checks only `order`, `ex/ey/ez` and `lx/ly/lz`. `order_min`, `dt`, `timemax`, `das_samples`, `das_gauge_length` and the SLS vector sizes (help says the sizes must match) are not checked. The string options `implem`, `method`, `mesh`, `anisotropy` and `das_type` are not validated here either.
- `SemProxyOptions`: `dt`, `timemax` and `taper_delta` are initialized from double literals (no `f` suffix), unlike the other float fields. The header has both an include guard and `#pragma once`.
- `SemProxyOptions::bind_cli`: `boundaries_size`, `srcx/y/z`, `rcvx/y/z` and the domain sizes are documented in meters, but nothing checks that the source and receiver lie inside the domain. The DAS azimuth and dip conventions (reference axes, sign) are not defined anywhere in this file.
- `src/main/fe/src/main.cc`, `main`: `ParseOptions` may call `exit()` on `--help` or invalid options after `Kokkos::initialize` (and MPI init) without finalizing either, and the `rank` and `size` values obtained from `InitMpi` are never used.
- `src/main/fe/src/main.cc`, `InitMpi`, `FinalizeMpi`, `ParseOptions`: non-static functions with external linkage in a driver file, and `ParseOptions` relies on `cxxopts` being included transitively through `sem_proxy_options.h`.
- `src/main/fe/src/main.cc`, `main`: `setenv` of `OMP_PROC_BIND` and `OMP_PLACES` unconditionally overrides values set by the user, and is called after MPI init, which may be too late if the OpenMP runtime is already initialized.
- `src/model/mesh/impl/builder/cartesian/include/cartesian_model_file_reader.h`, `CartesianModelFileReader`: the header uses `std::move` without including `<utility>`, and the closing `#endif` comment names a guard (`..._INCLUDE_MODEL_FILE_READER_H_`) different from the one defined (`..._INCLUDE_CARTESIAN_MODEL_FILE_READER_H_`).
- `CartesianModelFileReader::parse`: a section with count 0 leaves `count_` at 0, so the next section can set a different count and the "count mismatch" check is bypassed; the function also ends with a stray `;` after its body, and a `Model` header line found inside a value block is read as a non-numeric value.
- `CartesianModelFileReader::parse`: when the file ends right after a header, `count_line` is empty and the error reports it as a non-integer count. `std::stoul` accepts trailing garbage and negative numbers (`"-1"` wraps around), and `std::stod` accepts trailing garbage, so malformed values can pass unnoticed.
- `CartesianModelFileReader`: `parse` is private and called from the constructor, and the class is a public builder input but sits under `impl/`, so its file format and units are documented nowhere else; the units of each property are not stated by the reader.
- `src/model/mesh/impl/builder/cartesian/include/cartesian_params.h`: the header uses `std::string` without including `<string>`, and has no includes at all, so it relies on the includer.
- `CartesianParams`: the default constructor leaves `order`, `ex/ey/ez`, `lx/ly/lz`, `isModelOnNodes` and `isElastic` uninitialized, and the 9-argument constructor leaves `isAcoustoElastic`-related and global/origin fields at defaults, so global and origin fields must be set separately by the caller with no invariant tying them to the local sizes.
- `CartesianParams`: `DgSemBoundaryZ` and `acoustoElasticBoundaryZ` use inconsistent naming (PascalCase vs camelCase) and mix solver-specific coupling data into a generic mesh parameter struct.
- `src/model/mesh/impl/builder/cartesian/include/cartesian_partitioner.h`, `CartesianXPartitioner::partition`: when `size` exceeds `global.ex`, some ranks get `local.ex = 0` and an empty subdomain, with no check or error. `global.ex == 0` divides by zero in `dx`, and the `ScalarType` `rank * base_ex` product can overflow for large integer types.
- `CartesianXPartitioner::partition`: only the X-range is decomposed while the class is named as a generic Cartesian partitioner, and `<cmath>` is included but unused.
- `src/model/mesh/impl/builder/cartesian/include/cartesian_struct_boundary_classifier.h`, `CartesianStructBoundaryClassifier::classify`: calls unqualified `fabs` without a visible `<cmath>` qualification and applies it to `FloatType`, and `n_node` is not checked against `nx*ny*nz`, so an inconsistent count reads or writes out of range.
- `CartesianStructBoundaryClassifier`: `ScalarType` is only used to cast flag values that are then stored in a `vectorInt`, so a `ScalarType` other than `int` is silently converted back to `int`.
- `CartesianStructBoundaryClassifier::classify`: a node on a global z_max face with `free_surface_on_top` false but also on another global face gets Damping, and a Surface node on a lateral edge loses its Damping status, so the corner and edge rule is implicit.
- `src/model/mesh/impl/builder/cartesian/include/cartesian_struct_builder.h`, `CartesianStructBuilder::DgSemBoundaryZ_`: the constructor parameter is stored but read by nothing in the class; dead parameter.
- `CartesianStructBuilder` constructor: the member initializer list is not in declaration order (`ox_` and `global_o*_` are declared before `ex_`), which triggers `-Wreorder`; the 21-parameter positional constructor mixes local, global and coupling data and is easy to permute silently.
- `CartesianStructBuilder::getModel`, acousto-elastic branch: the fluid and solid vp/vs/rho values (1500/3400, 0/1963, 1020/2500) are hardcoded, the acoustic tensors are not driven by user data; if a model file is also given, its properties silently replace the two-layer values.
- `CartesianStructBuilder::getModel`, model file branch: only the vp extent is checked against the element count; rho and vs arrays are sized by the file count and never checked. A file lacking `Vp` gives a misleading "0 elements" error, and the debug `std::cout` of vp values reads a possibly device-resident view from the host.
- `CartesianStructBuilder::getModel`: `tol` divides by `ex_`, `ey_`, `ez_` unchecked (zero gives a division by zero), and `n_elem`, `n_node` are computed as `int` regardless of `ScalarType`.
- `src/model/mesh/impl/builder/cartesian/include/cartesian_struct_builder.h`: the header uses `std::string`, `std::shared_ptr`, `std::runtime_error`, `std::cout`, `std::to_string` and `std::move` without including `<string>`, `<memory>`, `<stdexcept>`, `<iostream>`, `<utility>`, relying on transitive includes.
- `src/model/mesh/impl/builder/cartesian/include/cartesian_unstruct_boundary_classifier.h`, `CartesianUnstructBoundaryClassifier::classify`: calls unqualified `fabs` on `FloatType`, and `n_node` is not checked against the size of the coordinate views. `ScalarType` is only used to cast flag values that are then stored in a `vectorInt`, so a `ScalarType` other than `int` is silently converted back to `int`.
- `CartesianUnstructBoundaryClassifier::classify`: the same corner and edge rule as the structured classifier is implicit: a Surface node on a lateral edge loses its Damping status. The coordinates are read on the host from views that may be device-resident.
- `src/model/mesh/impl/builder/cartesian/include/cartesian_unstruct_builder.h`, `CartesianUnstructBuilder::DgSemBoundaryZ_`: the constructor parameter is stored but read by nothing in the class; dead parameter.
- `CartesianUnstructBuilder` constructor: the member initializer list is not in declaration order (`ox_`, `global_o*_` are declared before `ex_`), which triggers `-Wreorder`. The default constructor leaves `ex_`, `ey_`, `ez_`, `lx_`..`lz_`, `order_` and the flags uninitialized, and the constructor is not `explicit`.
- `CartesianUnstructBuilder::initModels`: the acousto-elastic two-layer values differ between the node branch (rho 1020/2500) and the element branch (rho 1000/2000), so the same options give different models depending on `isModelOnNodes`. The element branch adds the origin `oz_` to find the layer, while the node branch uses `nodes_coords_z_`, which carries no origin (known zero offset), so the layer boundary is inconsistent across ranks. All values are hardcoded and not driven by user data.
- `CartesianUnstructBuilder::initModels`, model file branch: only the vp extent is checked against the element count, and the check comes after the debug `std::cout` that reads `model_vp_element_[n_element / 2]` and `[n_element - 1]` (out of range on a mismatch) from a possibly device-resident view. Rho, vs and other arrays are never checked. A file without `Vp` keeps the uniform vp, and file properties silently replace the two-layer acousto-elastic values.
- `CartesianUnstructBuilder::getCoordInOneDirection`, `initNodesCoords`: coordinates are computed and buffered as `float` although the class is templated on `FloatType` and `vectorReal` follows `real_t`; the `offset` argument is always 0, and the `global_i < nodes_x` guard is always true.
- `CartesianUnstructBuilder::getModel`: `tol` divides by `ex_`, `ey_`, `ez_` unchecked, and `n_node` is computed as `int` regardless of `ScalarType`. The header uses `std::string`, `std::shared_ptr`, `std::cout`, `std::to_string` and `std::min` without including their headers.
- `src/model/mesh/impl/builder/cartesian/src/cartesian_unstruct_builder.cc`, `CartesianUnstructBuilder<float, int>`: the only instantiation hardcodes `float`, so the builder cannot follow a `real_t = double` build.
- `src/model/mesh/impl/common/include/elasticity_utils.h`, `computeCTensor`: duplicates the VTI construction of `computeVTICoefficients` / `buildVTITensor` instead of reusing them, and mixes double literals (1.0, 2.0, 0.0) with `FloatType` so a `float` instantiation computes intermediate terms in double.
- `elasticity_utils.h`, `computeVTICoefficients`, `computeCTensor`: `sqrt` of `diff^2 + 2 vp^2 delta diff` is unchecked and yields NaN for strongly negative delta; `sqrt`, `cos` and `sin` are called unqualified without including `<cmath>`, relying on `data_type.h`.
- `elasticity_utils.h`, `buildIsotropicTensor`, `buildVTITensor`, `computeCTensor`: the `FloatType C[6][6]` parameters decay to pointers to `FloatType[6]`, so nothing checks the caller's array size, and the header has no include guard problem but no `<cmath>` / `common_macros.h` include of its own (`PROXY_HOST_DEVICE` comes via `data_type.h`).
- `elasticity_utils.h`, `computeIsotropicCoefficients`, `buildIsotropicTensor`, `computeVTICoefficients`, `buildVTITensor`: free function templates with generic names in the global namespace, in an `impl/common` directory.
- `src/model/mesh/impl/model/struct/include/gllpoints.h`, `GLLPoints::get`: an invalid order or index silently returns 0.0f, a valid-looking coordinate, instead of failing, so a wrong call goes undetected.
- `gllpoints.h`: the header is not self-contained: it uses `PROXY_HOST_DEVICE` without including `common_macros.h`. `<array>` and `<cstddef>` are included but unused, and `MAX_GLL_ORDER` is an unprefixed global constant.
- `GLLPoints::get`: the points are hardcoded as `float` literals (16 digits given) and do not follow `real_t`, so a `real_t = double` build gets single-precision nodes.
- `src/model/mesh/impl/model/struct/include/model_struct.h`, `ModelStruct::nodeCoord`: the branch `localIdx == Order` can never be true because `localIdx = nodeIdx % Order`; dead code. The `switch` on `dim` also has no default, and `dim` outside 0..2 indexes `nodeIdx` out of range.
- `ModelStructData`, `ModelStruct` constructor: `dx_/dy_/dz_` are copied into `lx_/ly_/lz_` (domain lengths) and the element size is then derived as `lx_ / ex_`. The `d` prefix suggests an element size, so the meaning is contradictory.
- `ModelStructData`: `ex_..ez_`, `dx_..dz_`, `isModelOnNodes_` and `isElastic_` have no default initializer. A default-constructed data object copied into a `ModelStruct` reads indeterminate values. `ModelStruct::ex_`, `ey_`, `ez_`, `ox_`.. are likewise uninitialized in the default constructor.
- `ModelStruct` constructor: the member initializer list is not in declaration order (`ox_` and `lx_` come before `isModelOnNodes_`, while the members are declared `ex_, nx_, lx_, hx_, ox_, ...`), which triggers `-Wreorder`. `boundaries_t_` and `face_connectivity_` are also initialized in a different order than declared.
- `ModelStruct::getCTensorOnElement`: reads `model_C_tensor_element_` without checking it was allocated. `initElasticityTensors` only allocates it for TTI, so calling it for Iso or VTI reads an empty view. It is also the only anisotropy path, and it ignores the per-element data (known hardcoded values).
- `ModelStruct::setQualityFactors`, `setModelNodeProps`: they write from host loops into `vectorReal` views that may be device-resident, and `setModelNodeProps` is not `PROXY_HOST_DEVICE` while the getters that read the same array are.
- `ModelStruct::getModelDelta*`, `getModelEpsilon*`, `getModelGamma*`: `ModelStructData` and `ModelStruct` have no storage for these, so they always return 0 and a heterogeneous anisotropic model cannot be expressed on structured meshes.
- `ModelStruct::getMinSpacing`: returns -1 for an order above 9 instead of failing; the tabulated factors are `float` literals and do not follow `FloatType`.
- `ModelStruct::elemOwner`, `elemNeighbor`, `localFaceOwner`, `localFaceNeighbor`: public helpers that are not overrides of `ModelApi`, so they are unreachable through the base interface (extends the duplicated-interface flag on face queries).
- `src/model/mesh/impl/model/unstruct/include/face_connectivity_unstruct.h`, `FaceConnectivityUnstruct::build`: the owner code `elem * 8 + lf` is cast to `int` while `elem` is `ScalarType`, so it overflows for more than about 268 million elements. `max_faces` is cast to `uint32_t` unchecked.
- `FaceConnectivityUnstruct::build`, Pass D: the neighbor dof buffer is a hardcoded `kMaxDofsPerFace = 100` with no check against `ndofs_per_face`, and a face dof with no match in the neighbor list leaves `face_perm` / `face_perm_inv` unset without any error.
- `FaceConnectivityUnstruct::build`: `face_id_of_bucket`, `elem_to_faces_` and the owner/neighbor tables are `vectorInt` / `arrayInt` while element and face ids are `ScalarType`, so a `ScalarType` wider than `int` is silently truncated. `FaceConnectivityUnstructData::n_faces` and the constructor argument are likewise `ScalarType` over `int` tables.
- `FaceConnectivityUnstruct(const FaceConnectivityUnstructData&)`: not `explicit`, so a data struct converts implicitly to a connectivity. Nothing checks that the injected table shapes agree with `n_faces` and `ndofs_per_face`.
- `FaceConnectivityUnstruct`, `FaceConnectivityUnstructData`: `FloatType` is unused by both classes.
- `FaceConnectivityUnstruct::getNumberOfFaces` and siblings: the constructor from data is marked `PROXY_HOST_DEVICE` although it copies Kokkos views, which is not meaningful on the device for a class holding a vtable (virtual base).
- `src/model/mesh/impl/model/unstruct/include/model_unstruct.h`, `ModelUnstructData`: the full constructor never initializes `ox_`, `oy_`, `oz_`, which `ModelUnstruct(const ModelUnstructData&)` copies, so the local origin is indeterminate. The unrelated `origin_x_/y_/z_` members (default 0) duplicate the notion of origin.
- `ModelUnstruct`: the default constructor leaves all scalar members uninitialized, and the constructor from data is not `explicit`. Its initializer list is not in declaration order (`n_points_per_element_` is declared before `isModelOnNodes_`, `model_phi_*` is initialized before `model_theta_*`), which triggers `-Wreorder`.
- `ModelUnstruct::faceNormal`: the `switch` on `CubicFace` has no default, so an invalid face leaves `n0`, `n1`, `n2` uninitialized. A degenerate face (norm below 1e-12) returns an unnormalized vector silently.
- `ModelUnstruct::getCTensorOnElement`: reads `model_C_tensor_element_` without checking it was allocated; `initElasticityTensors` only allocates it for TTI, so Iso and VTI read an empty view. `setModelNodeProps` and `setQualityFactors` write from host loops into views that may be device-resident.
- `ModelUnstruct`: the header uses `std::numeric_limits`, `std::runtime_error`, `fmin`, `sqrt` and `max` without including `<limits>`, `<stdexcept>`, `<cmath>`; `<array>` and `<map>` are included but unused. `getMaxSpeed` and `getMinSpacing` return a `-lowest`/`max` sentinel if no data is present in the element-only or node-only case, and `getMaxSpeed` is not host/device marked while `getMinSpacing` is.
- `src/model/mesh/impl/model/unstruct/src/model_unstruct.cc`, `ModelUnstruct` explicit instantiations: the four instantiations (`FloatType` float/double, `ScalarType` int/long) are compiled here, but the sibling builder (`CartesianUnstructBuilder<float, int>`) and the differentiators only instantiate `<float, int>`, so the `double` and `long` variants are never exercised by the rest of the code and are likely untested.
- `src/model/mesh/pywrap/include/bindings_builder.h`, `bind_cartesian_struct_builder` (second constructor): the global domain lengths are hardcoded to -1 and the global origin to 0 in the lambda, so Python cannot build a multi-rank subdomain. The -1 sentinel and the 20 positional arguments are easy to permute silently.
- `bind_cartesian_struct_builder`: the two constructors expose different subsets of the C++ constructor (model file and acousto-elastic parameters only in the second), and the acousto-elastic lambda returns `T` by value, which requires a movable builder.
- `bind_modelbuilderbase`, `bind_cartesian_struct_builder`, `bind_cartesian_unstruct_builder`: the header uses `std::shared_ptr` without including `<memory>`, and has both an include guard and `#pragma once`. Registration order (base before derived, params before unstruct builder) is an unchecked implicit requirement, as for the other bindings.
- `bind_cartesian_unstruct_params`: `CartesianParams` fields such as `origin_x/y/z` are exposed without the global sizes or the DG-SEM coupling fields, so a Python-built params object cannot describe a multi-rank subdomain either; the `.def(py::init<int, ...>)` takes `int` for `order` while the other counts follow `ScalarType`.
- `src/model/mesh/pywrap/include/bindings_face_connectivity.h`, `bindFaceConnectivityUnstruct`: the header uses `std::string` without including `<string>`, relying on transitive includes. The bindings of `FaceConnectivityUnstruct` are only registered for `ModelUnstruct<FloatType, ScalarType>` (via `build<MeshType>`), and the class must be registered after the model class is known to Python: an unchecked implicit ordering requirement.
- `bindFaceConnectivityUnstruct`: the data-struct property getters take `FaceConnData &` and return the Python view of the table, and the setters rebind the view without checking its shape against `n_faces` and `ndofs_per_face`, so a Python-built data object can be inconsistent.
- `src/model/mesh/pywrap/include/bindings_model.h`, `bind_modelstructdata`: only `ex`, `ey`, `ez`, `dx`, `dy`, `dz` are exposed, so Python cannot set `isModelOnNodes_`, `isElastic_` or the model arrays, which have no default initializer (known); a Python-built `ModelStructData` cannot describe a usable model.
- `bind_modelapi`: the binding exposes the getters (`get_model_delta_*`, `get_model_epsilon_*`, `get_model_gamma_*`, `get_model_theta_*`, `get_model_phi_*`) that are known to return 0 or truncated values, and `is_boundary_face` / `get_global_node_from_face` / `get_global_face` inherit the duplicated face interface: Python callers see the same inconsistencies with no warning.
- `bind_modelapi`: `set_quality_factors` and `init_elasticity_tensors` are bound on the base class although their effect depends on the implementation (device-resident views, allocation only for TTI); the constraints are not visible from Python.
- `src/model/mesh/pywrap/include/bindings_model.h`: the header has both an include guard and `#pragma once`, uses `std::shared_ptr` without including `<memory>`, and registration order (`ModelApi` before `ModelStruct`/`ModelUnstruct`, face connectivity before `ModelUnstructData`) is an unchecked implicit requirement. The binder templates are also unrelated to the `model` namespace only by convention, and `namespace py = pybind11;` is declared at global scope in a header.
- `bind_modelunstructdata`: the 30-argument constructor is bound with 22 positional `python_view_type_t<vectorReal>` views that are easy to permute silently. The Python-side constructor omits `ox`, `oy`, `oz`, so the origin is indeterminate (known for the C++ constructor).
- `src/model/mesh/pywrap/include/bindings_partionner.h`, `bind_cartesian_partitioner`: the header (and its include guard) is misspelled "partionner", and `namespace py = pybind11;` is declared at global scope in a header.
- `bind_cartesian_partitioner`: `partition` is bound without any return value policy or documentation of its C++ signature; the Python-side semantics of `rank` and `size` (and the empty-subdomain case when `size` exceeds the element count) are not visible from the binding.
- `src/model/mesh/pywrap/include/bindings_utils.h`: the header uses `std::runtime_error`, `std::invalid_argument` and `std::integral_constant` without including `<stdexcept>` and `<type_traits>`, relying on the includer. It has both an include guard and `#pragma once`.
- `order_suffix`, `orderDispatch`: the same unsupported-order error is reported with two different exception types (`std::runtime_error` versus `std::invalid_argument`). The supported range 1..9 is hardcoded in both and not tied to `MAX_GLL_ORDER`.
- `orderDispatch`: the return type is deduced from nine separate `return` statements, so the functor must return exactly the same type for every order; a functor returning an order-dependent type fails to compile with an unclear error. The name also duplicates the private `orderDispatch` helper of `differentiator_factory.cc`.
- `src/model/mesh/pywrap/include/bindings_utils.h`: the header uses `std::runtime_error`, `std::invalid_argument` and `std::integral_constant` without including `<stdexcept>` and `<type_traits>`, relying on the includer. It has both an include guard and `#pragma once`.
- `order_suffix`, `orderDispatch`: the same unsupported-order error is reported with two different exception types (`std::runtime_error` versus `std::invalid_argument`). The supported range 1..9 is hardcoded in both and not tied to `MAX_GLL_ORDER`.
- `orderDispatch`: the return type is deduced from nine separate `return` statements, so the functor must return exactly the same type for every order; a functor returning an order-dependent type fails to compile with an unclear error. The name also duplicates the private `orderDispatch` helper of `differentiator_factory.cc`.
- `src/model/mesh/pywrap/src/bindings.cpp`, `PYBIND11_MODULE(model, m)`: the C++ module is named `model` and its `__name__` is overwritten with the hardcoded string "pyfuntides.model", so the Python-visible name depends on the package layout and not on the build target.
- `src/model/mesh/pywrap/src/bindings.cpp`: the `orderDispatch` lambdas return `nullptr` only to satisfy the deduced return type, and each order loop instantiates 4 template combinations for 9 orders (36 classes per binder), which multiplies compile time.
- `src/model/mesh/pywrap/src/bindings.cpp`: `<cstdint>`, `<string>`, `<pybind11/numpy.h>` and `<pybind11/stl.h>` are included but nothing in this file uses them directly.
- `src/model/mesh/pywrap/src/bindings.cpp`: every binder is registered for `<double, *>` and `<*, long>`, while the builders and other back-ends are only instantiated for `<float, int>`; these Python classes may fail to link or are untested.
- `src/parallel/include/parallel_topology.h`: the header uses `size_t` without including `<cstddef>`, relying on transitive includes from `<map>` / `<vector>`.
- `src/parallel/include/parallel_topology.h`, `ParallelTopology`: nothing enforces that `myRank < numRanks` or that `sharedNodes` keys are valid ranks distinct from `myRank`; the struct is fully public with no invariant check.
- `src/solver/fe/DG-SEm/impl/acoustic/include/dg-sem_physics_traits_acoustic.h`, `DGSEMPhysicsTraits`: the struct is not templated on physics and its name does not say "acoustic", while a sibling `PhysicsTraits` in `solver::fe` also exists; no evidence in this file that anything consumes `kName`, `WavefieldType` or `RhsType`.
- `src/solver/fe/DG-SEm/impl/acoustic/include/dg-sem_rhs_acoustic.h`, `DGSEMRhsAcoustic`: the DG and SEM sources are built from the same `element` and `weights`, so both domains must have the same source element indices and weights. `getElement()` and `getWeights()` return only the DG member's data and, unlike `getTerm`, are not marked `override`, so they do not follow the per-component convention of the base interface.
- `DGSEMRhsAcoustic::getTerm`: any `i` other than 0 silently returns the SEM term, with no range check against `kNumRhsComponents`.
- `DGSEMRhsAcoustic`: a struct with public `m_`-prefixed data members that are also the way sub-solvers access their source, so the class is not encapsulated and the `m_` prefix contradicts its public status.
- `src/solver/fe/DG-SEm/impl/acoustic/include/dg-sem_wavefield_acoustic.h`, `DGSEMWavefieldAcoustic::getDGCurrentField`, `getSEMCurrentField`, `getDGPreviousField`, `getSEMPreviousField`: the index `i` is ignored and field 0 is always returned, so a wrong index goes undetected.
- `DGSEMWavefieldAcoustic`: a struct with public `m_`-prefixed data members, so the `m_` prefix contradicts the public status; `getNumFields`, `getFieldNames`, `swap` and `print` are not host/device marked and the struct has no base class in common with the sibling wavefield types.
- `DGSEMWavefieldAcoustic`: the constructor is not `explicit` (harmless with four arguments), and the four positional array arguments (two `arrayReal`, two `vectorReal`) are easy to permute silently.
- `src/solver/fe/DG-SEm/impl/common/include/dg-sem_solver.h`, `DGSEMsolver::getMassMatrixAcoustic`, `getDampingMatrix`, `getForceVector`: the exception messages say "not implemented for DG" although thrown by the coupled DG-SEM solver, unlike `getMassMatrixElastic` which says "DG-SEM coupling".
- `DGSEMsolver::m_penalty_factor_`: the SIPG penalty is duplicated here and in the DG sub-solver with no mechanism enforcing that the two stay equal (initialized to 12.0f in both places independently).
- `DGSEMsolver::m_mesh_`: `MESH_TYPE` is held by value, so the coupled solver keeps its own copy of the mesh in addition to the one passed to `computeFEInit`, and `m_face_connectivity_` is shared with the DG sub-solver only by convention.
- `dg-sem_solver.h`: the header uses `std::runtime_error` without including `<stdexcept>`, relying on transitive includes.
- `DGSEMsolver`: the destructor is declared `= default` in a class deriving from `Solver`; whether it is virtual depends on `Solver`, and the empty overrides (`computeForces`, `updateSolutionForward`, `updateSolutionBackward`) silently do nothing instead of failing if called.
- `src/solver/fe/DG-SEm/impl/common/include/dg-sem_solver_data.h`, `DGSEMsolverData`: a struct with public `m_`-prefixed members, so the prefix contradicts the public status; the constructor copies the wavefield and the source (Kokkos views, so shallow) with no check that they are consistent.
- `DGSEMsolverData::isDistributed`: a public flag with no documented meaning in this header and no visible reader.
- `src/solver/fe/DG-SEm/impl/common/include/dg-sem_solver_impl.h`, `DGSEMsolver::ApplyCoupling`: coordinates and local buffers (`faceCoords`, `dg_coords`, `sem_coords`, `normal_dg`, `stiff_dg_local`, `norm_dg`, `norm_sem`, `acc_*`) are hardcoded `float` while the integral callbacks and penalties use `real_t`, so the kernel cannot follow a `real_t = double` build (same float/real_t mismatch as the integral back-ends).
- `DGSEMsolver::ApplyCoupling`, `computeOneStep`: the coupled solver reaches into the DG and SEM sub-solvers through public members (`m_DG_solver_.m_stiff_local_`, `m_list_mode_`, `m_elem_list_`, `m_face_list_`, `m_n_*_list_`) and toggles `m_list_mode_` around the step; if a call between the two assignments throws, the DG sub-solver stays in list mode.
- `DGSEMsolver::TagElements`, `BuildDGInteriorFaceList`: `m_element_type_`, `DG_elem_list_`, `SEm_elem_list_`, `m_interface_face_indices_` and `SEm_node_list_` are filled by host loops through `operator[]` on views created by `allocateVector`, which may be device-resident.
- `DGSEMsolver::TagNodes`: the interface-face detection loop is duplicated verbatim (count pass and fill pass), and the `e >= nElem` guard in the kernel is redundant with the `RangePolicy` bound.
- `DGSEMsolver::computeFEInit`: `computeGlobalMassMatrixMasked` and `computeDampingMatrixMasked` are called on the SEM sub-solver after it has already assembled the full-mesh matrices, so the unmasked assembly is done and discarded.
- `DGSEMsolver::computeOneStep`: `dt` and `timeSample` are passed as `const float&` / `const int&` for scalars, and `dynamic_cast<DataType&>` throws `std::bad_cast` with no message on a mismatch.
- `src/solver/fe/DG-SEm/pywrap/include/bindings_dgsem_solver.h`, `bind_dgsem_rhs_acoustic`, `bind_dgsem_acoustic_data`: non-template, non-inline functions defined in a header, so including it in more than one translation unit violates the one-definition rule.
- `bind_dgsem_rhs_acoustic`, `bind_dgsem_acoustic_data`: registering them before their base classes (`Rhs`, `Solver::DataStruct`) fails at import time; this order is an unchecked implicit requirement. The header also has no include for `Kokkos::Experimental::python_view_type_t` and relies on the includer, and `DGSEMWavefieldAcoustic` is bound nowhere, so Python cannot build the wavefield argument of `DGSEMsolverData`.
- `bind_dgsem_rhs_acoustic`: the constructor takes four positional views (two `arrayReal`, one `vectorInt`, one `arrayReal`) that are easy to permute silently; the DG and SEM terms and `weights` have no documented shapes.
- `src/solver/fe/DG-SEm/pywrap/include/bindings_dgsem_wavefield.h`, `bind_dgsem_wavefield_acoustic`: non-template, non-inline function defined in a header, so including it in more than one translation unit violates the one-definition rule.
- `bind_dgsem_wavefield_acoustic`: the header uses `std::shared_ptr` without including `<memory>`, and registration order (the class must be bound before `DGSEMsolverData` uses it) is an unchecked implicit requirement. `namespace py = pybind11;` is declared at global scope in a header.
- `bind_dgsem_wavefield_acoustic`: the constructor takes four positional views (two `arrayReal`, two `vectorReal`) that are easy to permute silently, and the get_*_field bindings expose an index that the C++ getters ignore.
- `src/model/mesh/impl/model/struct/include/face_connectivity_struct.h`, `FaceConnectivityStruct::localFaceNeighbor`: the opposite face is obtained with `localFaceOwner(face_id) ^ 1`, which silently depends on the `CubicFace` enumerators being ordered as pairs (XMinus/XPlus, YMinus/YPlus, ZMinus/ZPlus); nothing checks it.
- `FaceConnectivityStruct`: the header has no include of its own for `PROXY_HOST_DEVICE`, and a default-constructed object (all counts 0) makes `getGlobalFace`, `getGlobalNodeFromFace` and `elemOwner` divide by zero. The constructor is marked `PROXY_HOST_DEVICE` although the class has a virtual base (same concern as `FaceConnectivityUnstruct`).
- `FaceConnectivityStruct`: `FloatType` is unused, and the members `ex_`..`nx_` are `ScalarType` while `order_` is `int`, so `order_ * ex_` mixes types.
- `src/parallel/include/topology_factory.h`, `TopologyFactory::createFromMesh`: the `catch (...)` around `getMinSpacing()` swallows every exception, and the tolerance is computed as `minDx * 1e-4` with a `double` literal, then stored in a `double` field. The `tol.auto_compute = false` assignment in the catch block has no effect because `tol` is a local copy that is not read afterwards.
- `TopologyFactory::createFromMesh`: `mesh.nodeCoord(i, 0)` returns the local coordinate, which for `ModelUnstruct` carries no subdomain origin (known), so comparing it to `origin_x` and `origin_x + domain_width_x` only works if the caller passes a consistent origin. The static `getMinSpacing` of unstructured meshes measures only element 0 (known), which feeds the auto tolerance.
- `TopologyFactory::createFromMesh`: the `topo.sharedNodes[...]` indices are cast with `static_cast<int>(i)` from `ScalarType`, so a wider `ScalarType` is silently truncated. The header has `using namespace utils;` at file scope, which pollutes every includer, and `<iostream>` is included but unused.
- `TopologyFactory`: a class with only one static template method and no state; the decomposition is hardcoded to the X direction (left/right neighbors are rank-1 and rank+1), consistent with `CartesianXPartitioner` but not stated by the name.
- `src/solver/fe/DG/impl/acoustic/include/dg_physics_traits_acoustic.h`, `DGPhysicsTraits`: the struct is not templated on physics and its name does not say "acoustic", while a sibling `PhysicsTraits` in `solver::fe` also exists; no evidence in this file that anything consumes `kName`, `WavefieldType` or `RhsType`.
- `src/solver/fe/DG/impl/acoustic/include/dg_wavefield_acoustic.h`, `DGWavefieldAcoustic::getCurrentField`, `getPreviousField`: the index `i` is ignored and the field is always returned, so a wrong index goes undetected.
- `DGWavefieldAcoustic`: a struct with public `m_`-prefixed data members, so the prefix contradicts the public status; `getNumFields`, `getFieldNames`, `swap` and `print` are not host/device marked, and the struct has no common base with the sibling wavefield types.
- `dg_wavefield_acoustic.h`: the header uses `std::swap`, `std::cout` and `std::endl` without including `<utility>` and `<iostream>`, relying on `data_type.h`.
- `src/solver/fe/DG/impl/common/include/dg_penalty.h`, `computeHexVolume`: the corner differences are scaled by 0.25 and the result by 8, so for the unit cube (corners at 0 and 1) the function returns 8 instead of 1. The factors look inconsistent with a volume of 8 * det at the center (a 0.125 factor would be expected for a [-1,1] parent element), which would make `h_f` and the SIPG penalty gamma off by a constant factor.
- `dg_penalty.h`: `computeFaceArea`, `computeHexVolume` and the penalty functions are non-template, non-inline function definitions in a header (`computeFaceArea`, `computeHexVolume`), so including it in more than one translation unit violates the one-definition rule; `sqrt` and `fabs` are used with float literals (0.5f, 8.0f) that do not follow `real_t = double`.
- `dg_penalty.h`, `computeSIPGPenaltyFromArea`: no check that `area` is nonzero, so a degenerate face gives a division by zero in `h_f`.
