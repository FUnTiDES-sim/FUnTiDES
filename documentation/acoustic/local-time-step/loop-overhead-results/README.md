# Surcoût de la boucle du pas de temps local (étape 1 du plan)

Mesures CPU du 9 novembre 2026 (Kokkos Serial, Release, cœur 0, `taskset -c 0`) sur le
cas gmsh à inclusion sur-raffinée de
[`gmsh-inclusion-stability-results`](../gmsh-inclusion-stability-results/README.md) :
7 répétitions en ordre aléatoire après 2 répétitions de chauffe, boucle de propagation
seule. Le scénario `Control_local_macro` (LTS sur le maillage témoin, sans élément fin)
mesure le surcoût de la boucle locale ; il a été ajouté au benchmark à cette occasion.

Sur le maillage raffiné (2 899 éléments), 144 éléments sont au petit pas, mais le support
de la correction compte **362 éléments** (12,5 %) : tous les éléments qui touchent un nœud
fin. Avec la même efficacité par élément, le gain est donc borné par 2 / 1,125 ≈ 1,78.

| Version | Surcoût sans zone fine | Boucle LTS | Global au petit pas | Gain | Erreur finale LTS |
|---|---:|---:|---:|---:|---:|
| Avant (`timing_summary_before.csv`) | +14,9 % | 6 247 ms | 9 494 ms | 1,52 | 4,41e-4 |
| a. Boucle allégée (`timing_summary_step1a.csv`) | +0,3 % | 5 340 ms | 9 468 ms | 1,77 | 4,01e-4 |
| b. + métrique précalculée (`timing_summary_step1b.csv`) | +0,6 % | **5 150 ms** | 9 401 ms | **1,83** | 4,07e-4 |

- **a.** Le masque de surface libre n'est appliqué que dans les éléments qui la touchent
  (repérés à l'initialisation) et `1/ρ` est calculé une fois par élément ;
  l'accélération fine `P M⁻¹ K u` n'est plus remise à zéro à chaque pas (seuls les nœuds
  fins sont écrits) ; la correction est soustraite directement dans `K u` par le noyau
  restreint, ce qui supprime un tampon, sa remise à zéro et une passe sur tous les nœuds.
- **b.** Les facteurs géométriques `w α B` des éléments du support sont calculés une fois
  par `configureLocalTimeStep` ; le noyau de correction ne recalcule plus jacobiennes et
  métriques. La raideur globale reste celle de `sem2d`, pour une comparaison équitable.

Les critères de l'étape 1 sont atteints : surcoût sans zone fine inférieur à 3 %, gain
d'au moins 1,8. La validation complète du cas passe (`test_result_step1b.json`) ; seules
les erreurs du calcul local changent, par l'ordre des additions en simple précision :
erreur finale 4,07e-4 au lieu de 4,41e-4, récepteurs 4,15e-4 au lieu de 4,21e-4, bilan
d'énergie identique (5,47e-5). Les 70 tests unitaires LTS / SEM 2D / gmsh passent sur CPU,
112 sur GPU (CUDA).

Reproduction, depuis la racine du dépôt (le benchmark écrit ses CSV dans le répertoire
du cas : utiliser une copie contenant `bilayer.msh`, `materials.csv`,
`timing_parameters.txt` et le `snapshots.csv` de la validation) :

```bash
cmake --build build-lts-cpu --target sem2d_gmsh_inclusion_timing
env -u LD_LIBRARY_PATH taskset -c 0 \
  build-lts-cpu/tests/validation/sem-local-time-step-2d/sem2d_gmsh_inclusion_timing \
  documentation/acoustic/local-time-step/gmsh-structured-results <copie> <copie>/timing_parameters.txt 7
```
