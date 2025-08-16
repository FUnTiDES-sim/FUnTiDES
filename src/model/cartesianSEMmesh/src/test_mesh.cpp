#include <cartesianSEMmesh.hpp>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

// Old implementation of mesh
void globalNodesList(int size, int order, vector<vector<int>> &nodesList) {
  int nx = size * order + 1;
  int nz = nx;
  for (int k = 0; k < size; k++) {
    for (int j = 0; j < size; j++) {
      for (int i = 0; i < size; i++) {
        int n0 = i + j * size + k * size * size;
        int offset = i * order + j * order * nx + k * order * nx * nz;
        for (int m = 0; m < order + 1; m++) {
          for (int n = 0; n < order + 1; n++) {
            for (int l = 0; l < order + 1; l++) {
              int dofLocal =
                  l + n * (order + 1) + m * (order + 1) * (order + 1);
              int dofGlobal = offset + l + n * nx + m * nx * nz;
              nodesList[n0][dofLocal] = dofGlobal;
            }
          }
        }
      }
    }
  }
}

std::vector<float> getCoordInOneDirection(const int &order, const int &nCoord,
                                          const int &h, const int &nElement) {
  std::vector<float> coord(nCoord);
  std::vector<float> xi(order + 1);
  switch (order) {
  case 1:
    xi[0] = -1.;
    xi[1] = 1.;
    break;
  case 2:
    xi[0] = -1.;
    xi[1] = 0.;
    xi[2] = 1.;
    break;
  case 3:
    static constexpr float sqrt5 = 2.2360679774997897;
    xi[0] = -1.0;
    xi[1] = -1. / sqrt5;
    xi[2] = 1. / sqrt5;
    xi[3] = 1.;
    break;
  case 4:
    static constexpr float sqrt3_7 = 0.6546536707079771;
    xi[0] = -1.0;
    xi[1] = -sqrt3_7;
    xi[2] = 0.0;
    xi[3] = sqrt3_7;
    xi[4] = 1.0;
    break;
  case 5:
    static constexpr float sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
    static constexpr float sqrt__7_mins_2sqrt7__ = 1.30709501485960033;
    static constexpr float sqrt_inv21 = 0.218217890235992381;
    xi[0] = -1.0;
    xi[1] = -sqrt_inv21 * sqrt__7_plus_2sqrt7__;
    xi[2] = -sqrt_inv21 * sqrt__7_mins_2sqrt7__;
    xi[3] = sqrt_inv21 * sqrt__7_mins_2sqrt7__;
    xi[4] = sqrt_inv21 * sqrt__7_plus_2sqrt7__;
    xi[5] = 1.0;
    break;
  default:
    break;
  }
  int i = nElement;
  float x0 = i * h;
  float x1 = (i + 1) * h;
  float b = (x1 + x0) / 2.;
  float a = b - x0;
  for (int j = 0; j < order + 1; j++) {
    coord[j] = a * xi[j] + b;
  }
  return coord;
}

void nodesCoordinates(vector<vector<float>> &nodeCoordsX,
                      vector<vector<float>> &nodeCoordsZ,
                      vector<vector<float>> &nodeCoordsY, int size, int order,
                      float elemSize) {
  std::vector<float> coordX(order + 1);
  std::vector<float> coordY(order + 1);
  std::vector<float> coordZ(order + 1);

  for (int n = 0; n < size; n++) {
    coordZ = getCoordInOneDirection(order, order + 1, elemSize, n);
    for (int m = 0; m < size; m++) {
      coordY = getCoordInOneDirection(order, order + 1, elemSize, m);
      for (int l = 0; l < size; l++) {
        coordX = getCoordInOneDirection(order, order + 1, elemSize, l);
        int e = l + m * size + n * size * size;
        for (int k = 0; k < order + 1; k++) {
          for (int j = 0; j < order + 1; j++) {
            for (int i = 0; i < order + 1; i++) {
              nodeCoordsX[e][i + (order + 1) * j +
                             k * (order + 1) * (order + 1)] = coordX[i];
              nodeCoordsZ[e][i + (order + 1) * j +
                             k * (order + 1) * (order + 1)] = coordZ[k];
              nodeCoordsY[e][i + (order + 1) * j +
                             k * (order + 1) * (order + 1)] = coordY[j];
            }
          }
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  const int size = 3;
  const int order = 5;
  const float domainSize = 200.;
  const float elemSize = domainSize / size;
  const int numberOfElement = size * size * size;
  const int numberOfNodes =
      (size * order + 1) * (size * order + 1) * (size * order + 1);
  // init new mesh
  CartesianParams<float, int> params{ order, size, size, size, domainSize, domainSize, domainSize };
  CartesianSEMmesh<float, int, order> mesh(params);

  // security checks
  assert(numberOfElement == mesh.getNumberOfElements());
  assert(numberOfNodes == mesh.getNumberOfNodes());

  // 1) Checking if Element To node map is valide
  cout << "Testing if element map is right, testing first corner...\n";
  vector<vector<int>> oldGlobalNodeList(numberOfElement, vector<int>(numberOfNodes, 0));
  globalNodesList(size, order, oldGlobalNodeList);
  for (int z = 0; z < size; z++) {
    for (int y = 0; y < size; y++) {
      for (int x = 0; x < size; x++) {
        int element = x + y * size + z * size * size;
        assert(element < size * size * size);
        for (int nz = 0; nz < order + 1; nz++) {
          for (int ny = 0; ny < order + 1; ny++) {
            for (int nx = 0; nx < order + 1; nx++) {
              int nodeId = nx + ny * (order+1) + nz * (order + 1) * (order+1);
              if (oldGlobalNodeList[element][nodeId] != mesh.globalNodeIndex(element, nx, ny, nz))
              {
                cout << "Element (" << element << ") (" << x << ", " << y << ", " << z
                     << ") at node (" << nx << ", " << ny << ", " << nz << ") give ";
                cout << mesh.globalNodeIndex(element, nx, ny, nz) << " against " << oldGlobalNodeList[element][nodeId]
                     << endl;
                // return 1;
              }
            }
          }
        }
      }
    }
  }
  cout << "... SUCCESS.\n";

  // 2) Checking Coord
  cout << "Testing node coordinate against old Mesh api...\n";
  vector<vector<float>> oldCoodX(numberOfElement,
                                 vector<float>(numberOfNodes, 0));
  vector<vector<float>> oldCoodY(numberOfElement,
                                 vector<float>(numberOfNodes, 0));
  vector<vector<float>> oldCoodZ(numberOfElement,
                                 vector<float>(numberOfNodes, 0));
  nodesCoordinates(oldCoodX, oldCoodZ, oldCoodY, size, order, elemSize);

  // Imprecision tolerance for float
 const double EPSILON = 2;

  // Checking coordinates
  for (int z = 0; z < size; z++) {
    for (int y = 0; y < size; y++) {
      for (int x = 0; x < size; x++) {
        int elemIndex = x + y * size + z * size * size;
        for (int nz = 0; nz < order + 1; nz++) {
          for (int ny = 0; ny < order + 1; ny++) {
            for (int nx = 0; nx < order + 1; nx++) {
              int localNodeIndex = nx + ny * (order + 1) + (order + 1) * (order + 1) * nz;
              int globalNodeIndex = mesh.globalNodeIndex(elemIndex, nx, ny, nz);

              float oldCoordxN = oldCoodX[elemIndex][localNodeIndex];
              float oldCoordyN = oldCoodY[elemIndex][localNodeIndex];
              float oldCoordzN = oldCoodZ[elemIndex][localNodeIndex];

              float newCoordxN = mesh.nodeCoord(globalNodeIndex, 0);
              float newCoordyN = mesh.nodeCoord(globalNodeIndex, 1);
              float newCoordzN = mesh.nodeCoord(globalNodeIndex, 2);

              if (fabs(oldCoordxN - newCoordxN) > EPSILON ||
                  fabs(oldCoordyN - newCoordyN) > EPSILON ||
                  fabs(oldCoordzN - newCoordzN) > EPSILON) {
                std::cout << "Mismatch at element " << elemIndex << " node " << localNodeIndex << std::endl;
                std::cout << "Old: (" << oldCoordxN << ", " << oldCoordyN << ", " << oldCoordzN << ") "
                          << "New: (" << newCoordxN << ", " << newCoordyN << ", " << newCoordzN << ")" << std::endl;
                return 0;
              }
            }
          }
        }
      }
    }
  }

  std::cout << "All coordinates match within EPSILON." << std::endl;
  cout << "...SUCCESS on X axis.\n";

  return 0;
}
