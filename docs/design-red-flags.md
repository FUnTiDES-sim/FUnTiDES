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
