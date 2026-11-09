# Interface Gmsh → FUnTiDES

Les lecteurs produisent un `ModelUnstructData<float, int>` directement utilisable
par `ModelUnstruct<float, int>`. Le lecteur 3D alimente les solveurs acoustiques
SEM 3D existants ; le solveur `sem-local-time-step-2d` reste en 2D.

## Maillages acceptés

- `readGmsh2D` : quadrilatères Q4 dans le plan XZ, groupes physiques de surfaces
  pour les matériaux et de lignes pour les frontières.
- `readGmsh3D` : hexaèdres Hex8 (type Gmsh 5), groupes physiques de volumes pour
  les matériaux acoustiques `(vp, rho)` et de surfaces Q4 pour les frontières.
- Format natif : **MSH 2.2 ASCII**. Les identifiants peuvent être non contigus.
- Géométrie linéaire Q4 / trilinéaire Hex8 ; les nœuds GLL sont reconstruits à
  l'ordre SEM demandé, indépendamment de l'ordre géométrique de Gmsh.
- Maillage conforme : chaque face extérieure doit avoir un groupe physique de
  frontière ; chaque face intérieure appartient à deux cellules et ne doit pas
  être exportée comme frontière physique. Les lignes et points 3D sont ignorés.

Les maillages mixtes, les éléments géométriques d'ordre supérieur et les
jacobiennes non positives sont rejetés. Le contrôle 3D échantillonne la jacobienne
sur une grille de référence de `(2*p+1)^3` points : il ne certifie pas sa positivité
partout dans une cellule très déformée. Les cellules inversées doivent être
corrigées dans le maillage source. La priorité des nœuds de surface libre est
conservée aux intersections avec les frontières absorbantes.

Le partage des nœuds GLL utilise la topologie des sommets, arêtes et faces,
y compris lorsque deux cellules ont des orientations locales différentes.
Il ne fusionne pas des points sur la seule base de leur proximité géométrique.

## Exemple reproductible

Depuis la racine de ce dépôt, avec un Python disposant du module `gmsh` :

```bash
python3 src/model/mesh/impl/builder/gmsh/examples/generate_hex_box.py /tmp/box.msh
cmake --build build-lts-cpu --target gmsh_sem3d_inspect
build-lts-cpu/src/model/mesh/impl/builder/gmsh/gmsh_sem3d_inspect /tmp/box.msh /tmp/box.csv 2
```

Le générateur crée 8 hexaèdres, un matériau de tag 10 et deux frontières :
1 = surface libre supérieure, 2 = absorption sur les cinq autres côtés.
À l'ordre 2, le lecteur construit 125 nœuds GLL, dont 25 de surface libre
et 73 absorbants. Le CSV a les colonnes `physical_tag,vp,rho`.

Le convertisseur existant détecte maintenant la dimension et accepte les Hex8 :

```bash
python3 src/solver/impl/sem-local-time-step-2d/validation/convert_gmsh_mesh.py input41.msh output22.msh
```

L'option `--plane xy` est réservée au 2D. La conversion conserve les groupes
physiques exportés par Gmsh ; elle ne transforme pas des tétraèdres en hexaèdres.

## Utilisation C++

```cpp
#include "gmsh_3d_reader.h"

auto data = model::readGmsh3D("box.msh", 2,
    {{10, {1500.f, 1000.f}}},
    {{1, model::BoundaryFlag::Surface}, {2, model::BoundaryFlag::Damping}});
auto mesh = std::make_shared<model::ModelUnstruct<float, int>>(data);
mesh->buildFaceConnectivity();
// Passer *mesh à computeFEInit du solveur SEM 3D acoustique,
// configuré à l'ordre 2 et avec les propriétés par élément.
```

Le lecteur est disponible sans `COMPILE_SEM_2D`. L'import initialise un matériau
acoustique par élément ; les propriétés élastiques/anisotropes et le
partitionnement MPI de fichiers Gmsh ne sont pas ajoutés par cette interface.

## Utilisation Python

Avec les bindings compilés (`ENABLE_PYWRAP=ON`) :

```python
from pyfuntides import model

data = model.read_gmsh_3d(
    "box.msh", 2,
    {10: model.GmshAcousticMaterial(1500.0, 1000.0)},
    {1: model.BoundaryFlag.Surface, 2: model.BoundaryFlag.Damping},
)
mesh = model.ModelUnstruct_f32_i32(data)
mesh.build_face_connectivity()
```

`read_gmsh_2d` expose également le lecteur 2D existant.

## Vérification

```bash
cmake --build build-lts-cpu --target gmsh_3d_reader_test gmsh_2d_reader_test
build-lts-cpu/tests/units/model/gmsh_3d_reader_test
build-lts-cpu/tests/units/model/gmsh_2d_reader_test
```

Les tests 3D vérifient le partage d'une face tournée aux ordres 1, 2, 3, 4 et 9,
les groupes physiques, les erreurs de format/topologie/géométrie et, lorsque
`COMPILE_SEM` est actif, les masses et forces d'un champ non constant comparées
au constructeur cartésien ainsi que la conservation d'un champ constant après
une mise à jour temporelle.
