# Validation visuelle : trois couches de vitesses différentes

Les figures et les CSV proviennent du calcul C++ réalisé par
[`runThreeLayerValidation()`](three_layer_case.h). Le même calcul est exécuté
par le test `SEM2DLTS.ThreeLayersAdaptToVelocityAndMatchSmallGlobalStep`.
Python lit les CSV et réalise uniquement les tracés.

## Configuration

- Domaine : 6 000 × 3 000 m ; `z` vers le haut.
- SEM acoustique d’ordre 2, éléments carrés de 100 m dans les trois couches :
  1 800 éléments et 7 381 nœuds. Densité uniforme égale à 1.
- Impulsion gaussienne initiale centrée en `(3000, 2500)` m, écart-type 160 m,
  dérivée temporelle initiale nulle. Aucune source pendant le calcul.
- Bords rigides, sans amortissement. Les réflexions sur les bords font donc
  partie de ce cas de validation.
- Durée : 1,8 s ; stabilisation `nu=0.01`, CFL de partition 0,4.

| Couche | Intervalle z | Vp | Pas estimé | Pas retenu |
|---|---|---|---|---|
| Supérieure | 2 000–3 000 m | 1 500 m/s | 13,33 ms | 12 ms |
| Intermédiaire | 1 000–2 000 m | 2 250 m/s | 8,89 ms | 6 ms |
| Inférieure | 0–1 000 m | 3 000 m/s | 6,67 ms | 6 ms |

Le maillage est identique dans les trois couches : le changement de niveau
temporel est donc lié à la vitesse. La version actuelle possède deux niveaux,
avec un rapport 2. Les deux couches rapides partagent le petit pas ; il n’y a
pas encore trois pas distincts. Les nœuds partagés avec une zone à petit pas
sont inclus dans cette zone. Le choix par estimation géométrique ne remplace
pas la vérification de stabilité.

## Résultats mesurés sur CPU

La référence utilise le même maillage et un pas global de 1,5 ms. Le calcul
local effectue 150 grands pas, la référence 1 200 petits pas.

- Erreur finale relative L², pondérée par la masse : **1,62003 %**.
- Dérive relative maximale de l’énergie discrète : **2,6993 × 10⁻⁷**.
- 600 éléments à grand pas et 1 200 éléments à petit pas.

Ces résultats comparent deux intégrations temporelles sur la même discrétisation
spatiale. Ils ne mesurent pas l’erreur spatiale face à une solution exacte et
ne constituent pas une mesure de gain de performance.

## Comparaison des temps de calcul

Mesures CPU sur Intel Core i9-13900HX, Kokkos Serial, compilation Release GCC
12.3, un thread fixé au processeur logique 0. Chaque scénario est échauffé deux
fois, puis mesuré 11 fois dans un ordre mélangé. La boucle est encadrée par des
`Kokkos::fence()`. Elle exclut initialisation, diagnostics d’énergie, sorties,
échantillonnage et copies vers le host. La construction du maillage est commune
et exclue ; la préparation du solveur et des champs est mesurée séparément.

| Intégration | Nombre de pas | Boucle médiane | Intervalle interquartile | Préparation médiane | Erreur relative |
|---|---:|---:|---:|---:|---:|
| Local 12 / 6 ms | 150 grands pas | 65,41 ms | 65,37–65,48 ms | 1,98 ms | 1,643 % |
| Global 9 ms | 200 | 46,17 ms | 46,16–46,22 ms | 1,28 ms | 1,133 % |
| Global 6 ms | 300 | 74,48 ms | 74,43–74,51 ms | 1,28 ms | 0,499 % |
| Global 3 ms | 600 | 149,10 ms | 148,98–149,23 ms | 1,28 ms | 0,119 % |
| Global 1,5 ms | 1 200 | 290,95 ms | 290,81–291,10 ms | 1,29 ms | 0,0247 % |

Les erreurs de ce tableau utilisent une référence plus fine à **0,75 ms**, sur
le même maillage. Elles diffèrent donc légèrement des valeurs de la validation
visuelle, dont la référence est à 1,5 ms.

Sur ce cas, le prototype local réduit le temps de boucle de **12,2 %** par
rapport au global à 6 ms, mais avec une erreur plus élevée. Le global à 9 ms
est à la fois plus rapide et plus précis que le prototype local. Il n’y a donc
pas de gain global de performance/précision démontré pour le LTS sur ce modèle.
Deux tiers des éléments utilisent le petit pas, ce qui limite le potentiel du
LTS ; le prototype effectue aussi des corrections et des parcours de vecteurs
supplémentaires. Ces mesures ne préjugent pas des performances sur GPU ni sur
des maillages où la zone à petit pas serait plus petite.

L’[extension du 7 novembre à huit modèles](fraction-results/README.md) fait varier
la proportion fine de 0 à 100 % et affine les pas pour un budget d’erreur finale
de 2 %. Sur les huit modèles, les candidats locaux retenus restent plus lents
que les candidats globaux. Les données, la figure et les limites de cette
comparaison sont conservées dans `fraction-results/`.

Les mesures brutes sont dans `results/timing_samples.csv`, les médianes dans
`results/timing_summary.csv`, le protocole et l’environnement dans
`results/timing_protocol.json` et `results/timing_environment.json`.

## Figures

1. [Modèle et pas locaux](results/01_modele_et_pas_locaux.png)
2. [Propagation à quatre instants](results/02_propagation_trois_couches.png)
3. [Comparaison à la référence et champ de différence](results/03_comparaison_reference.png)
4. [Erreur et énergie en fonction du temps](results/04_erreur_et_energie.png)
5. [Temps de calcul et précision](results/05_temps_calcul.png)

Les amplitudes de pression sont en unités arbitraires. Les cartes de propagation
partagent la même échelle de couleur. Dans la comparaison finale, l’écart possède
une échelle distincte, indiquée sur sa barre de couleur.

## Reproduire

Depuis la racine du worktree, après la configuration CPU décrite dans le
[README du composant](../README.md) :

```bash
cmake --build build-lts-cpu --target sem2d_lts_three_layers_validation \
    sem_solver_2d_local_time_step_test -j 4
env -u LD_LIBRARY_PATH build-lts-cpu/tests/units/solver/sem_solver_2d_local_time_step_test \
    --gtest_filter=SEM2DLTS.ThreeLayersAdaptToVelocityAndMatchSmallGlobalStep
env -u LD_LIBRARY_PATH build-lts-cpu/src/solver/impl/sem-local-time-step-2d/sem2d_lts_three_layers_validation \
    src/solver/impl/sem-local-time-step-2d/validation/results
../FUnTiDES-install/tpl-venv/bin/python \
    src/solver/impl/sem-local-time-step-2d/validation/plot_three_layers.py
```

Le script de tracé nécessite NumPy et Matplotlib. Les fichiers `elements.csv`,
`fields.csv`, `history.csv` et `metadata.json`, dans `results/`, conservent les
valeurs ayant servi aux figures. Les commandes remplacent les résultats de ce cas.

Pour reconstruire le benchmark :

```bash
cmake --build build-lts-cpu --target sem2d_lts_three_layers_benchmark -j 4
```

Pour le relancer avec le même cœur et tracer les résultats :

```bash
env -u LD_LIBRARY_PATH taskset -c 0 \
    build-lts-cpu/src/solver/impl/sem-local-time-step-2d/sem2d_lts_three_layers_benchmark \
    src/solver/impl/sem-local-time-step-2d/validation/results 11
../FUnTiDES-install/tpl-venv/bin/python \
    src/solver/impl/sem-local-time-step-2d/validation/plot_three_layer_timings.py
```

Adapter le cœur si le processeur logique 0 n’est pas accessible dans la session.
Les valeurs du tableau ci-dessus décrivent la mesure enregistrée et ne sont pas
mises à jour automatiquement par une nouvelle exécution.
