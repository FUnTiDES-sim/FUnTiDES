# `src/solver` — layout and type composition

## Layout

```
solver/
├── api/include/   solver.h, wavefield.h, rhs.h, physics_traits.h — the pure contract
├── impl/
│   ├── common/                 shared across every method
│   │   ├── include/, src/      solver_factory.{h,cc}
│   │   ├── acoustic/include/   rhs_acoustic.h, physics_info_acoustic.h
│   │   ├── elastic/include/    rhs_elastic.h, physics_info_elastic.h
│   │   └── acoustoelastic/include/  rhs_acoustoelastic.h
│   ├── sem/         common/ + acoustic/ + elastic/ + acoustoelastic/
│   ├── dg/          common/ + acoustic/
│   ├── dgsem/       common/ + acoustic/
│   └── dg_padaptive/ common/ + acoustic/
└── pywrap/          the single `solver` Python module (like gradient/pywrap)
```

| Tree | CMake option | Content |
|---|---|---|
| `api` | (always) | Pure interfaces common to every method: `solver.h` (base class), `wavefield.h`, `rhs.h`, `physics_traits.h` |
| `impl/common` | (always) | Construction entry point `solver_factory.{h,cc}` and the per-physics types shared across methods (`rhs_*.h`, `physics_info_*.h`) |
| `impl/sem` | `COMPILE_SEM` | SEM solvers: `common/` skeleton plus one subtree per physics (`{include,src}`) |
| `impl/dg` | `COMPILE_DG` | Acoustic DG solver (base of the other variants): `common/` + `acoustic/` |
| `impl/dgsem` | `COMPILE_DG_SEM` (implies `COMPILE_DG`) | Coupled DG-SEM (DG + SEM domains), composite reusing `impl/sem/` + `impl/dg/` |
| `impl/dg_padaptive` | `COMPILE_DG_PADAPTIVE` (implies `COMPILE_DG`) | p-adaptive DG (pMin/pMax subdomains), composite reusing `impl/dg/` |
| `pywrap` | `ENABLE_PYWRAP` | The pyfuntides `solver` module; per-method bindings guarded by the `COMPILE_*` definitions |

Every method compiles into the single `fun_solver` library. All
`<method>/{common,<physics>}/include` roots are on its include path, which
supports the cross-tree includes documented below. The factory and the
bindings guard their per-method parts with the `COMPILE_*` definitions, so
a build with `COMPILE_SEM=OFF` compiles the shared boundary only.

## Wavefield / RHS type composition

Data-type deduplication goes through **composition**, not inheritance:
each variant reuses the base tree's types and defines only its own
accessors. No data type is copied.

```
WavefieldAcoustic        (impl/sem/acoustic, nodal SEM: 2 vectorReal views)
DGWavefieldAcoustic      (impl/dg, per-element DG: 2 arrayReal views)
 ├── DGSEMWavefieldAcoustic    (dgsem) = DGWavefieldAcoustic + WavefieldAcoustic
 └── DGPAdaptiveWavefieldAcoustic (dg_padaptive) = DGWavefieldAcoustic x2 (pMin/pMax)

RhsAcoustic              (impl/common/acoustic, SEM: term/element/weights)
 ├── DGSEMRhsAcoustic    (dgsem) = 2 x RhsAcoustic (DG and SEM domains, shared element)
 └── DGPAdaptiveRhsAcoustic (dg_padaptive) = 2 x RhsAcoustic (pMin/pMax, distinct weights)
```

The base DG wavefield/RHS is reused as-is by the variants
(`DGSEMWavefieldAcoustic::getDGCurrentField()` delegates to
`DGWavefieldAcoustic::getCurrentField(0)`).

## Physics traits and the physics descriptor

Two description layers serve different consumers:

- **`PhysicsTraits`** (`api/include/physics_traits.h`, specialized per method
  tree, e.g. `impl/sem/acoustic/physics_traits_acoustic.h`) is the *type bundle*:
  `kName`, `WavefieldType`, `RhsType`. It binds concrete storage types and is
  consumed by a method's solver header. Each variant keeps its own traits
  file (`DGPhysicsTraits`, `DGSEMPhysicsTraits`, `DGPAdaptivePhysicsTraits`):
  their include sets differ (each variant exists only under its
  `COMPILE_DG_*` option, a merge would force `#ifdef` guards) and their
  struct names are meaningful to their solver headers. A new physics adds
  one struct per variant in the variant's own file.
- **`PhysicsInfo`** (`impl/common/include/physics_info.h`, specializations in
  `impl/common/<physics>/physics_info_<physics>.h`) is the
  *method-agnostic descriptor*: only `kName`, `kNumFields`, `kNumRhs` — the
  component counts intrinsic to the physics. Methods that need nothing
  beyond these counts (the DG variants) consume `PhysicsInfo` instead of a
  method's `PhysicsTraits`, so they do not depend on SEM storage types.

## pywrap bindings

The `solver` Python module is built once from `pywrap/`; per-method binding
headers stay separate because their signatures differ structurally (DG
`arrayReal` views vs SEM `vectorReal` views, per-domain getters). Merging
them would change the exposed Python API for a marginal line count.

## Factory

`solver_factory::createSolver` is the single construction entry point: it
dispatches (method, implementation, mesh, physics, order) to `SEMsolver`,
`SEMsolverAcoustoElastic`, `DGsolver`, `DGSEMsolver` or `DGPAdaptiveSolver`
according to the compile options. Its declaration header
(`impl/common/include/solver_factory.h`) is the public boundary common to
every method — it depends only on `solver.h` and the enums; the definition
(`impl/common/src/solver_factory.cc`) guards each method with the
`COMPILE_*` definitions.

## The `sem_solver_impl.h` fragments

The out-of-line definitions of `SEMsolver` live in thematic headers, all
re-included in one fixed order by the umbrella `sem_solver_impl.h`:

- `sem_solver_step.h` — `computeFEInit` and the one-step drivers
  (`computeForces`, `updateSolution*`, `reset*`, `applyRHSTerm`, the
  `computeElementContributions` dispatcher);
- `sem_solver_kernels.h` — the `detail::` capability probes, the z-deformed
  geometry helpers and the attenuation dispatcher, shared by both physics;
- `impl/sem/acoustic/sem_solver_acoustic_impl.h` — acoustic element kernels
  (GEMM/flat/team, z-deformed variants, acoustic attenuation);
- `impl/sem/elastic/sem_solver_elastic_impl.h` — elastic element kernels
  (Iso/Vti/Tti, flat/team/z-deformed variants, elastic attenuation, TTI
  tensor precomputation: `computeCMatrix`, `computeTtiCompact`);
- `sem_solver_update.h` — `updateFieldsForward/Backward`;
- `sem_solver_init.h` — setup and initialisation: global mass matrix,
  damping, array allocation, sponge values, diagnostics;
- `sem_solver_element_subset.h` — entry points on element/node subsets:
  masked assemblies, per-node-list contributions and field updates.

Include the fragment you need when adding code; the umbrella is what the
explicit-instantiation units include. The fragment boundary is also where a
future geometry pre-computation (today recomputed every step by the face
kernels) can plug in without touching the stepping.

## Per-physics hierarchy

Each physics lives in its method subtree — `impl/sem/acoustic/` and
`impl/sem/elastic/`, each with its types (wavefield, rhs, traits), its
implementation fragment (`sem_solver_{acoustic,elastic}_impl.h`) and its
explicit instantiations in `src/`. Besides the kernels above they carry the
named types:

- `AcousticSolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>` — a
  thin `final` subclass of `SEMsolver<..., kAcoustic>`;
- `ElasticSolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>` — the
  same for `kElastic` (Iso/Vti/Tti anisotropy stays a runtime model choice,
  as in the underlying template).

These names exist for readability; the factory remains the only
construction path. Dispatch stays compile-time (template parameter) — no
virtual method is added on a GPU path.

## The acousto-elastic fork

`SEMsolverAcoustoElastic` (`impl/sem/acoustoelastic/`) is **not** a mergeable
duplicate of `SEMsolver`: it derives from `Solver` directly (it does not
instantiate `SEMsolver`) and its specificity is physical, not structural —
element/node classification (`TagElements`, `TagNodes`), interface
properties (`SetInterfaceNodeProps`, `SaveInterfaceUnm1`) and the
bidirectional acoustic↔elastic coupling (`ApplyInterfaceCoupling`,
`ApplyCouplingAcousticToElastic`, `ApplyCouplingElasticToAcoustic`). These
steps have no equivalent in `SEMsolver`; factoring the shared mechanics
(mass, damping, update) would rework `SEMsolver` kernels for a ~150-line
gain at the highest physical risk in the repository. Decision: the fork
stays, and this paragraph is the rationale.
