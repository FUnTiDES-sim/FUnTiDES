# `src/solver` — layout et composition des types

## Arbres

| Arbre | Option CMake | Contenu |
|---|---|---|
| `SEM/api` + `SEM/impl` + `SEM/pywrap` | (toujours) | Solveurs SEM (acoustique, élastique, acoustoélastique) ; cible `fun_solver`, fabrique `solver_factory::createSolver` |
| `DG` | `COMPILE_DG` | Solveur DG acoustique (base des autres variantes) |
| `DG-SEm` | `COMPILE_DG_SEM` (implique `COMPILE_DG`) | DG-SEM couplé (domaines DG + SEM) |
| `DG_padaptive` | `COMPILE_DG_PADAPTIVE` (implique `COMPILE_DG`) | DG p-adaptatif (sous-domaines pMin/pMax) |

Toutes les variantes compilent dans la même bibliothèque `fun_solver` ; les
répertoires `impl/*/include` sont tous sur le chemin d'inclusion du target, ce
qui autorise les inclusions inter-arbres documentées ci-dessous.

Historique : la couche `fe/` (`src/solver/fe/...`) existait pour cohabiter
avec les solveurs `fd`, aujourd'hui supprimés ; elle a été retirée sur revue
(2026-10-19) : l'arbre SEM vit dans `SEM/`, les variantes DG à la racine, et
le namespace `solver::fe` est devenu `solver` (convention dossier =
namespace, comme `src/gradient`).

Historique : les répertoires `DG/api` et `DG-SEm/api` d'origine (un
`CMakeLists.txt` orphelin sans `add_subdirectory` correspondant, sans
`include/`) étaient morts et ont été supprimés. Réintroduits en
`refactor/dg-boundaries`, cette fois peuplés : chaque variante DG porte
désormais son `api/include/` avec sa face publique — l'en-tête de classe
(consommé par `SEM/impl/common/src/solver_factory.cc`), l'en-tête de données
(consommé par `src/main/fe`), et `dg_penalty.h` pour `DG` — déclaré comme
répertoire d'inclusion PUBLIC de `fun_solver` et installé par
`install_funtides_headers`, symétrique de `SEM/api`. Les `*_impl.h` et les
en-têtes par physique restent dans `impl/` (consommateurs à feu blanche :
tests unitaires). Les inclusions se font par nom nu : aucun consommateur
n'a changé. Le point d'entrée générique reste `SEM/api/include/solver.h`.

## Composition des types wavefield / RHS

La déduplication des types de données passe par la **composition**, pas par
l'héritage : chaque variante réutilise les types de l'arbre de base et ne
définit que ses propres accès. Aucun type de données n'est copié.

```
WavefieldAcoustic        (SEM/impl/acoustic, SEM nodal : 2 vues vectorReal)
DGWavefieldAcoustic      (DG, DG par élément : 2 vues arrayReal)
 ├── DGSEMWavefieldAcoustic   (DG-SEm) = DGWavefieldAcoustic + WavefieldAcoustic
 └── DGPAdaptiveWavefieldAcoustic (DG_padaptive) = DGWavefieldAcoustic + DGWavefieldAcoustic (pMin/pMax)

RhsAcoustic              (SEM/impl/acoustic, SEM : term/element/weights)
 ├── DGSEMRhsAcoustic    (DG-SEm) = 2 x RhsAcoustic (domaines DG et SEM, élément partagé)
 └── DGPAdaptiveRhsAcoustic (DG_padaptive) = 2 x RhsAcoustic (pMin/pMax, weights distincts)
```

Le wavefield/RHS DG de base est réutilisé tel quel par les variantes
(`getDGCurrentField()` de `DGSEMWavefieldAcoustic` délègue à
`DGWavefieldAcoustic::getCurrentField(0)`).

## Traits par variante

Chaque variante possède un *type bundle* de trois alias
(`kName`, `WavefieldType`, `RhsType`) consommé par son en-tête de solveur :

- `SEM/impl/acoustic` : spécialisation `PhysicsTraits<physicType::kAcoustic>` du
  gabarit déclaré dans `SEM/api/include/physics_traits.h` (dispatch par physique) ;
- `DG` : `DGPhysicsTraits` (`dg_physics_traits_acoustic.h`) ;
- `DG-SEm` : `DGSEMPhysicsTraits` (`dg-sem_physics_traits_acoustic.h`) ;
- `DG_padaptive` : `DGPAdaptivePhysicsTraits` (`dg_padaptive_physics_traits_acoustic.h`).

Ces fichiers sont volontairement restés distincts : contrairement aux
wavefields/RHS, leurs jeux d'inclusions diffèrent (chaque variante n'existe que
sous son option `COMPILE_DG_*`, une fusion imposerait des gardes `#ifdef`), et
leurs noms de structure sont porteurs (consommés par les en-têtes de solveur
respectifs). Une nouvelle physique ajoute une structure par variante dans le
fichier de sa variante.

## Bindings pywrap

Les bindings des variantes (`*/pywrap/include/bindings_*`) restent séparés
par variante : leurs signatures diffèrent structurellement (vues `arrayReal` DG
vs `vectorReal` SEM, getters spécifiques par domaine), et les fusionner
changerait l'API Python exposée pour un gain de lignes marginal.

## Fabrique

`solver_factory::createSolver` est le seul point d'entrée : il dispatch (méthode, implémentation, maillage, physique, ordre) vers `SEMsolver`, `SEMsolverAcoustoElastic`, `DGsolver`, `DGSEMsolver`
ou `DGPAdaptiveSolver` selon les options de compilation. Son en-tête de
déclaration, `solver_factory.h`, vit dans `SEM/api/include` (frontière publique,
ne dépend que de `solver.h` et des énumérations) ; l'implémentation reste dans
`SEM/impl/common/src/solver_factory.cc`.

## Découpage de `sem_solver_impl.h`

Les définitions out-of-line de `SEMsolver` sont réparties en en-têtes
thématiques, tous ré-inclus dans l'ordre original par le parapluie
`sem_solver_impl.h` (aucun TU ne change) :

- `sem_solver_step.h` — `computeFEInit` et les drivers de pas (`computeForces`,
  `updateSolution*`, `reset*`, `applyRHSTerm`, dispatcheur
  `computeElementContributions`) ;
- `sem_solver_kernels.h` (dans `impl/common/`) — les sondes de capacité
  `detail::`, les helpers de géométrie z-deformée, le dispatcheur
  d'atténuation, partagés par les deux physiques, et les constructeurs de
  tenseur TTI (`computeCMatrix`, `computeTtiCompact`) ;
- `impl/physics/acoustic/sem_solver_acoustic_kernels.h` — noyaux d'éléments
  acoustiques (GEMM/flat/team, variantes z-deformées, atténuation acoustique) ;
- `impl/physics/elastic/sem_solver_elastic_kernels.h` — noyaux d'éléments
  élastiques (Iso/Vti/Tti, variantes flat/team/z-deformées, atténuation
  élastique, précalcul des tenseurs TTI) ;
- `sem_solver_update.h` — `updateFieldsForward/Backward` ;
- `sem_solver_init.h` — setup et initialisation : matrice de masse globale,
  amortissement, allocation des tableaux, sponge, diagnostics (renommé depuis
  `sem_solver_assembly.h`, dont le nom promettait de l'assemblage) ;
- `sem_solver_element_subset.h` — points d'entrée sur sous-ensemble
  d'éléments/nœuds : assemblages masqués, contributions et mises à jour de
  champs par liste de nœuds.

Le découpage est purement mécanique (chaque extraction vérifiée être une
partition exacte de l'original, réassemblage octet-identique). C'est aussi la
couture où un futur pré-calcul de géométrie (aujourd'hui recalculée chaque pas
par les noyaux de faces) pourra se brancher sans toucher au stepping.

## Hiérarchie par physique (`impl/physics/`)

Depuis la migration en politiques (voir `documentation/physics-policy-design.md`),
chaque physique possède son répertoire : `impl/physics/acoustic/` et
`impl/physics/elastic/`, chacun avec sa racine d'inclusion ajoutée à
`fun_solver` (convention des arbres par variante). Outre les noyaux ci-dessus,
ils portent les types nommés :

- `AcousticSolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>` —
  sous-classe `final` mince de `SEMsolver<..., kAcoustic>` ;
- `ElasticSolver<ORDER, INTEGRAL_TYPE, MESH_TYPE, IS_MODEL_ON_NODES>` —
  idem pour `kElastic` (l'anisotropie Iso/Vti/Tti reste un choix d'exécution
  du modèle, comme pour le template sous-jacent).

Ces noms existent pour la lisibilité ; la fabrique reste le seul chemin de
construction. Le dispatch reste compile-time (paramètre template) — aucune
méthode virtuelle n'est ajoutée sur un chemin GPU.

## Fourche acoustoélastique

`SEMsolverAcoustoElastic` (`SEM/impl/acoustoelastic/`) n'est **pas** une
duplication fusionnable de `SEMsolver` : elle hérite directement de `Solver`
(elle n'instancie pas `SEMsolver`) et sa spécificité est physique, pas
structurelle — classification des éléments/nœuds (`TagElements`, `TagNodes`),
propriétés d'interface (`SetInterfaceNodeProps`, `SaveInterfaceUnm1`) et
application du couplage bidirectionnel acoustique↔élastique
(`ApplyInterfaceCoupling`, `ApplyCouplingAcousticToElastic`,
`ApplyCouplingElasticToAcoustic`). Ces étapes n'ont pas d'équivalent dans
`SEMsolver` ; factoriser la mécanique commune (masse, amortissement, update)
reviendrait à remanier les noyaux de `SEMsolver` pour un gain de ~150 lignes,
au risque physique le plus élevé du dépôt. Décision : la fourche reste, et le
présent paragraphe en est la motivation.
