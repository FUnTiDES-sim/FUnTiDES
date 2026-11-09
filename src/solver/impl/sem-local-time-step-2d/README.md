# SEM acoustique 2D avec pas de temps local

Développement sur la branche **`sem-local-timesetp-2d`**, dans le worktree
`FUnTiDES-sem-local-timesetp-2d/`. Le répertoire principal `FUnTiDES/` reste sur
`dev-2d` pour les développements de Claude. Ce composant est placé à côté du
solveur de référence [`sem2d/`](../sem2d/).

Le [document mathématique et d’implémentation](../../../../documentation/acoustic/local-time-step/documentation-mathematique-et-implementation.md)
présente les équations, la preuve de conservation de l’énergie, la correspondance
avec le code et l’architecture envisagée pour des pas indépendants non multiples.
Il contient les instantanés et les sismogrammes du
[test à petites inclusions](../../../../documentation/acoustic/local-time-step/inclusion-results/README.md).
L’[extension avec excitation en surface et quatre bords absorbants](../../../../documentation/acoustic/local-time-step/surface-absorbing-results/README-quatre-bords-absorbants.md)
ajoute la comparaison au témoin rigide, les traces différentielles de l’inclusion
et la vérification de la dissipation d’énergie.
La [configuration courante à surface libre](../../../../documentation/acoustic/local-time-step/free-surface-results/README.md)
remplace le bord supérieur absorbant par une pression nulle, avec excitation
sous la surface ; elle donne les temps de calcul global et local.

## Diagnostics et partition à deux niveaux

[`buildSEM2DLocalTimeStepPlan()`](include/sem_2d_local_time_step_plan.h) prépare
une partition fixe avec un grand pas `coarse_dt` et un petit pas `coarse_dt / 2`.
Il est compilé avec `COMPILE_SEM_2D=ON` et intégré à `fun_solver`.

Pour chaque élément, un noyau Kokkos lit les données sur le device et calcule :

- le plus petit espacement entre nœuds GLL adjacents, dans les deux directions ;
- la vitesse maximale de l’élément, ou de ses nœuds si le modèle est nodal ;
- une estimation `cfl * min_node_spacing / max_vp`.

Cette estimation sert au diagnostic et à la sélection des zones. **Elle ne
constitue pas une borne spectrale de stabilité**, notamment pour des éléments
distordus ou des matériaux hétérogènes. La stabilité de l’intégrateur local
devra être établie et testée séparément.

Les éléments dont l’estimation est inférieure au grand pas sont classés dans la
zone à petit pas. Si même `coarse_dt / 2` dépasse une estimation, la construction
échoue : il faut réduire le grand pas, ou ajouter ultérieurement d’autres niveaux.

Les listes produites sont triées par indice du maillage :

| Champ | Contenu |
|---|---|
| `element_steps` | Espacement minimal, vitesse maximale et estimation du pas par élément |
| `coarse_elements`, `fine_elements` | Éléments des deux zones |
| `fine_nodes` | Tous les nœuds des éléments à petit pas, y compris les nœuds partagés |
| `coarse_nodes` | Nœuds restants |
| `interface_nodes` | Nœuds appartenant à la fois à un élément à petit pas et à un élément à grand pas |
| `fine_support_elements` | Tous les éléments incidents à un nœud de `fine_nodes`, y compris les voisins de la zone à grand pas |

La dernière liste permet de préparer l’application de l’opérateur assemblé à un
champ porté par les nœuds de la zone à petit pas. La limiter aux seuls
`fine_elements` omettrait des contributions des voisins.

Exemple C++, après initialisation de Kokkos et construction d’un
`model::ModelUnstruct<float, int>` de dimension 2 :

```cpp
#include "sem_2d_local_time_step_plan.h"

auto plan = solver::buildSEM2DLocalTimeStepPlan(mesh, /*cfl=*/0.3f, /*coarse_dt=*/0.001f);
// Inspecter plan.element_steps et la proportion de plan.fine_elements.
// Aucune mise à jour du champ n’est effectuée par cette fonction.
```

Reconstruire la partition après un changement de géométrie, de vitesses ou de
pas de temps. Comme pour les solveurs, les modifications faites sur le host
doivent être synchronisées vers le device avant cet appel.

## Intégrateur local à deux niveaux

[`SEMsolver2DLocalTimeStep<ORDER, IS_MODEL_ON_NODES>`](include/sem_solver_2d_local_time_step.h)
implémente maintenant un premier intégrateur LF-LTS stabilisé, d’ordre 2 en
temps, pour un problème **sans source temporelle forcée**. Les conditions
absorbantes du premier ordre sont prises en charge par un amortissement
diagonal centré en temps. Il réutilise la masse
et la mise à jour de Verlet du SEM 2D, avec un opérateur de force modifié.

La formule est la spécialisation au rapport 2 du schéma de
[Grote, Michel et Sauter, section 2.4](https://arxiv.org/html/2005.13350v4#S2.SS4).
Avec `M` la masse diagonale, `K` la raideur et `P` la projection sur les nœuds
à petit pas, elle s’écrit :

```text
delta = 1 + nu / 4
a = (2 delta² - 1) / (16 delta²)
B u = K u - a ΔT² K P M⁻¹ K u
u_next = 2 u - u_prev - ΔT² M⁻¹ B u
```

Chaque grand pas nécessite une application de `K` sur tout le domaine et une
application restreinte à `fine_support_elements`. Les contributions des éléments
voisins sont donc présentes aux nœuds partagés. Les étapes intermédiaires sont
évaluées par cette correction polynomiale ; elles ne constituent pas des champs
directement disponibles pour des sorties à `ΔT/2`.

Le paramètre de stabilisation `nu` vaut `0.01` par défaut et doit appartenir à
`[0, 0.5]`. `nu=0` permet de retrouver le schéma non stabilisé. Le grand pas est
fixé lors de la configuration et doit rester identique dans les appels suivants.
Le plan géométrique ne garantit toujours pas la condition spectrale de stabilité.

### Utilisation C++

Initialiser Kokkos et construire au préalable un maillage acoustique 2D d’ordre 2,
avec propriétés par élément. Ses bords peuvent être rigides (condition naturelle,
`InteriorNode`), à pression nulle (`Surface`) ou absorbants du premier ordre
(`Damping`). Les maillages du constructeur
cartésien habituel portent par défaut des bords absorbants : ils sont maintenant
acceptés. Leur matrice de dissipation doit être finie et non négative.

```cpp
#include "sem_solver_2d_local_time_step.h"

solver::SEMsolver2DLocalTimeStep<2, false> lts;
lts.computeFEInit(mesh, {0.f, 0.f, 0.f}, false, 0.f);
lts.configureLocalTimeStep(dt, /*cfl=*/0.3f, /*nu=*/0.01f);

int const n = mesh.getNumberOfNodes();
vectorReal previous("previous", n), current("current", n);
// Remplir current avec le champ initial avant initializePreviousField().
solver::WavefieldAcoustic wavefield(previous, current);
solver::RhsAcoustic rhs(arrayReal("no source", 0, 0), vectorInt("no source element", 0),
                        arrayReal("no source weights", 0, 0));
solver::SEMsolverDataAcoustic data(wavefield, rhs);
lts.initializePreviousField(data);  // vitesse initiale nulle ; ne pas échanger les champs ici
double const initial_energy = lts.discreteEnergy(data);
for (int step = 0; step < number_of_steps; ++step) {
  lts.computeOneStep(dt, step, data);
  data.swapWavefields();  // un échange par grand pas
}
double const final_energy = lts.discreteEnergy(data);
```

`initializePreviousField(data, velocity)` accepte aussi une vitesse temporelle
initiale nodale. Cette vitesse est la dérivée temporelle de la pression, distincte
de la vitesse de propagation `vp`. Le démarrage utilise le même opérateur modifié
que la récurrence. L’énergie discrète est évaluée entre les deux niveaux stockés :

```text
E = 1/2 [(u-u_prev)ᵀ M (u-u_prev) / ΔT² + u_prevᵀ B u]
```

Sans amortissement, elle est conservée à l’arrondi près. Avec une matrice
diagonale de dissipation `D`, elle décroît de
`dt * v_center^T D v_center` par pas, où
`v_center = (u_next - u_prev) / (2*dt)`. Sa positivité dépend de la condition de
stabilité réelle. Les nœuds à pression nulle sont projetés à zéro dans les opérateurs
et les champs. Les tailles et les recouvrements mémoire des champs sont vérifiés.

`applyStiffness(field, result)` applique l’opérateur spatial à un champ quelconque.
`applyFineSupportStiffness(field, result)` sélectionne le voisinage local et
`applyEffectiveStiffness(field, result)` applique `B`. Chaque appel remplace son
résultat et utilise, pour le dernier, des tableaux de travail internes.

Le solveur est instancié une seule fois par couple ordre/emplacement du modèle,
dans des fichiers générés par CMake, pour éviter la duplication des lambdas CUDA.
Il est utilisable directement en C++, mais n’est pas encore exposé dans la fabrique
des solveurs ni dans les interfaces Python.

## Tests

Le binaire `sem_2d_local_time_step_plan_test` vérifie les cas homogènes, le
raffinement local, les variations de vitesse, les modèles nodaux, les mises à
jour du modèle sur le device, les voisinages et les entrées invalides. Les
diagnostics sont testés aux ordres 1, 2, 4 et 9, indépendamment de l’ordre maximal
compilé pour les intégrateurs.

Le binaire `sem_solver_2d_local_time_step_test` vérifie l’opérateur spatial face
au SEM de référence, les contributions des voisins y compris aux sommets, la
réduction au Verlet global sans zone à petit pas, une récurrence de Chebyshev
indépendante, les surfaces libres, la conservation de l’énergie et la convergence
temporelle d’un pulse traversant les deux frontières temporelles. Un cas compare
aussi la stabilité du LTS à la divergence de Verlet global au même grand pas.

Dans une compilation CPU configurée avec `COMPILE_SEM_2D=ON` :

```bash
cmake --build build-lts-cpu --target sem_2d_local_time_step_plan_test \
    sem_solver_2d_local_time_step_test sem_solver_2d_test -j 4
env -u LD_LIBRARY_PATH ctest --test-dir build-lts-cpu \
    -R 'SEM2DLocalTimeStepPlan|SEM2DLTS|SemSolver2DTest' --output-on-failure
```

La configuration CPU de développement est séparée de l’installation CUDA
existante. Aucun fichier de cette installation n’est remplacé.

Pour reproduire cette configuration locale depuis la racine du worktree. Kokkos
Serial est installé dans `../FUnTiDES-install-lts-cpu`, hors de `/tmp` (effacé au
redémarrage) et à part de l’installation CUDA :

```bash
cmake -S ../FUnTiDES-TPL/external/kokkos -B build-kokkos-serial -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/../FUnTiDES-install-lts-cpu" \
    -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF -DKokkos_ENABLE_CUDA=OFF \
    -DKokkos_ENABLE_TESTS=OFF -DKokkos_ENABLE_EXAMPLES=OFF
cmake --build build-kokkos-serial -j 4
cmake --install build-kokkos-serial
cmake -S . -B build-lts-cpu -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
    -DKokkos_DIR="$PWD/../FUnTiDES-install-lts-cpu/lib/cmake/Kokkos" \
    -Dbenchmark_DIR="$PWD/../FUnTiDES-install/lib/cmake/benchmark" \
    -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST="$PWD/../FUnTiDES/build-2d/_deps/googletest-src" \
    -DCOMPILE_SEM_2D=ON -DCOMPILE_DG=OFF -DCOMPILE_DG_SEM=OFF -DCOMPILE_DG_PADAPTIVE=OFF \
    -DENABLE_PYWRAP=OFF -DFUNTIDES_ENABLE_ADIOS2=OFF \
    -DMAX_SOLVER_ORDER=2 -DMAX_DIFFERENTIATOR_ORDER=1
```

Cette configuration réutilise les sources GoogleTest déjà présentes localement
pour éviter un téléchargement.

Les programmes et scripts de validation sont dans
[`tests/validation/sem-local-time-step-2d/`](../../../../tests/validation/sem-local-time-step-2d/)
(cibles CMake à la demande, non compilées par défaut). Les rapports, figures et résumés
sont dans [`documentation/acoustic/local-time-step/`](../../../../documentation/acoustic/local-time-step/),
avec la note mathématique. Les gros fichiers de sortie (sismogrammes `receivers.csv`,
instantanés `snapshots.csv`, champs et maillages exportés) ne sont pas versionnés
(`documentation/acoustic/local-time-step/.gitignore`) : les programmes les régénèrent.

### Résultats de validation CPU

**86 tests réussis** : 27 tests de partition, 43 tests de l’intégrateur local et
16 tests du SEM 2D de référence. Les intégrateurs sont testés aux ordres 1 et 2,
avec propriétés nodales et élémentaires. Huit tests vérifient aussi l’équivalence
au global amorti, le démarrage avec vitesse initiale et le bilan de dissipation.
Les cas ci-dessous emploient des unités
normalisées avec une vitesse de propagation de 1.

| Cas | Résultat observé |
|---|---|
| Énergie sur 1 500 grands pas, deux zones | Dérive relative maximale inférieure à `1.8e-6` sur les quatre configurations |
| Pulse traversant les frontières temporelles, durée 6 | Erreurs relatives pondérées par la masse : `3.7685e-3`, `9.21863e-4`, `2.22434e-4` pour `ΔT=0.12`, `0.06`, `0.03` |
| Même partition aux trois résolutions temporelles | Rapports d’erreur `4.09` et `4.14`, compatibles avec l’ordre 2 |
| Cas à `ΔT=0.25` avec raffinement local | Énergie conservée après 500 pas locaux ; amplitude de Verlet global supérieure à `1e5` après 50 pas au même grand pas |

La référence du pulse utilise le même maillage et Verlet avec un pas global de
`0.0075`. Le paramètre de partition `cfl` est réduit avec `ΔT` pour conserver la
même zone locale dans l’étude de convergence. Ces tests ne constituent pas une
preuve de stabilité pour tout maillage ou tout choix de pas.

Le diagnostic et les instanciations de l’intégrateur d’ordre 1 par élément et
d’ordre 2 nodal compilent aussi avec CUDA 12.9 pour `sm_89`. L’exécution GPU et
les performances restent à vérifier, le pilote n’étant pas accessible dans la
session. Les tailles de champs temporaires restent globales dans ce prototype.

### Cas à trois couches et figures

Le [cas de validation à trois couches](../../../../documentation/acoustic/local-time-step/validation-trois-couches.md) utilise un maillage
uniforme et des vitesses de 1 500, 2 250 et 3 000 m/s. Les pas retenus sont
respectivement de 12, 6 et 6 ms : les deux couches rapides partagent le petit pas
de la version à deux niveaux. À 1,8 s, l’erreur relative face au calcul global à
1,5 ms est de 1,62 %, avec une dérive relative maximale d’énergie de `2.7e-7`.
Les figures, CSV et commandes de reproduction sont dans `documentation/acoustic/local-time-step/`.

L’[étude de précision et de performances du 7 novembre](../../../../documentation/acoustic/local-time-step/fraction-results/README.md)
couvre huit modèles avec 0 à 100 % d’éléments au petit pas. Après raffinement
des pas pour une erreur finale plafonnée à 2 %, aucun candidat local retenu
n’est plus rapide que le candidat global du même modèle. Le témoin sans zone
fine présente déjà un surcoût d’environ 14 % pour la boucle du prototype.

## Plan d’amélioration

Le cas gmsh à inclusion sur-raffinée montre un gain : **9,6 s contre 15,1 s**
(1,57 fois) face au calcul global au petit pas, seul pas stable sur ce maillage,
avec 5 % d’éléments au petit pas. Le LTS est rentable quand le raffinement est
imposé par la géométrie (petite inclusion, couche mince, topographie) ; sur un
maillage adapté à la longueur d’onde, le rapport 2 n’apporte rien (étude du
7 novembre). Le gain idéal du cas à inclusion est d’environ 2 / 1,07 ≈ 1,9 ;
le surcoût de la boucle (14 % sans zone fine) explique l’essentiel de l’écart.

0. **Remise en état** — *fait le 9 novembre* : Kokkos Serial hors de `/tmp`,
   sorties volumineuses non versionnées, commit, rebase sur `dev-2d` (calcul
   de l’éponge indépendant du pas).
1. **Surcoût de la boucle** : profiler séparément la raideur globale, la raideur
   restreinte, les remises à zéro et la correction ; restreindre remises à zéro
   et correction aux nœuds fins et à leur voisinage (listes précalculées dans
   le plan) ; fusionner la projection `M⁻¹` et la correction dans le noyau
   restreint. Critères : surcoût sans zone fine inférieur à 3 %, gain sur
   l’inclusion gmsh d’au moins 1,8, erreurs inchangées.
2. **Domaine de rentabilité, CPU et GPU** : exécuter tests et benchmarks sur
   GPU (occupation du noyau restreint quand la zone fine est petite) ; mesurer
   le gain sur des maillages gmsh à raffinement géométrique en fonction de la
   proportion fine (1 à 30 %), le global prenant le pas imposé par sa stabilité.
   Résultat attendu : une règle d’emploi.
3. **Partition fondée sur une borne spectrale** : estimer `λmax(Mₑ⁻¹ Kₑ)` par
   élément (itération de puissance sur la matrice élémentaire) pour classer les
   éléments sans `cfl` arbitraire ; vérification globale optionnelle de
   `M⁻¹ B` par itération de puissance. Critère : partition sûre sur les maillages
   distordus, sans calibration spectrale séparée comme pour l’inclusion.
4. **Sources et récepteurs** : terme source dans le schéma, évalué aux temps
   intermédiaires dans la zone fine ; échantillonnage des récepteurs aux temps
   physiques, avec interpolation dans la zone fine. Critère : un tir Ricker
   identique à `sem2d` à la tolérance de l’étude de convergence.
5. **Intégration** : fabrique des solveurs, interfaces Python, lecteur gmsh,
   éponge ; option de pas local dans les tirs de makutu ; cas applicatif avec
   couche lente mince ou topographie maillée finement.
6. **Plusieurs niveaux** (rapports 4, 8, …, section 10 de la documentation
   mathématique) : multiplie le gain quand le raffinement géométrique est fort.
   Critère : trois pas distincts dans le cas à trois couches, énergie et ordre 2
   conservés, gain accru sur l’inclusion.
7. **Plus tard** : adjoint pour la FWI (le schéma est auto-adjoint : `B`
   symétrique), MPI compatible avec la partition, 3-D.

Les sources temporelles forcées, l’éponge, les champs adjoints et les données
distribuées provoquent actuellement une erreur explicite. Le solveur SEM 2D
existant et sa boucle temporelle sont conservés.
