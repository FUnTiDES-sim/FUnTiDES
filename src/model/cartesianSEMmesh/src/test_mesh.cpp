#include <cartesianSEMmesh.hpp>
#include <cassert>
#include <iostream>
#include <vector>

using namespace std;

// Old implementation of mesh
void globalNodesList(int size, int order, vector<vector<int>> &nodesList) {
  int nx = size * order + 1;
  int nz = nx;
  for (int j = 0; j < size; j++) {
    for (int k = 0; k < size; k++) {
      for (int i = 0; i < size; i++) {
        int n0 = i + k * size + j * size * size;
        int offset = i * order + k * order * nx + j * order * nx * nz;
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
    static constexpr double sqrt5 = 2.2360679774997897;
    xi[0] = -1.0;
    xi[1] = -1. / sqrt5;
    xi[2] = 1. / sqrt5;
    xi[3] = 1.;
    break;
  case 4:
    static constexpr double sqrt3_7 = 0.6546536707079771;
    xi[0] = -1.0;
    xi[1] = -sqrt3_7;
    xi[2] = 0.0;
    xi[3] = sqrt3_7;
    xi[4] = 1.0;
    break;
  case 5:
    static constexpr double sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
    static constexpr double sqrt__7_mins_2sqrt7__ = 1.30709501485960033;
    static constexpr double sqrt_inv21 = 0.218217890235992381;
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
    coordY = getCoordInOneDirection(order, order + 1, elemSize, n);
    for (int m = 0; m < size; m++) {
      coordZ = getCoordInOneDirection(order, order + 1, elemSize, m);
      for (int l = 0; l < size; l++) {
        coordX = getCoordInOneDirection(order, order + 1, elemSize, l);
        int e = l + m * size + n * size * size;
        for (int k = 0; k < order + 1; k++) {
          for (int j = 0; j < order + 1; j++) {
            for (int i = 0; i < order + 1; i++) {
              nodeCoordsX[e][i + (order + 1) * j +
                             k * (order + 1) * (order + 1)] = coordX[i];
              nodeCoordsZ[e][i + (order + 1) * j +
                             k * (order + 1) * (order + 1)] = coordZ[j];
              nodeCoordsY[e][i + (order + 1) * j +
                             k * (order + 1) * (order + 1)] = coordY[k];
            }
          }
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  const int size = 5;
  const int order = 2;
  const float domainSize = 200.;
  const float elemSize = domainSize / size;
  const int numberOfElement = size * size * size;
  const int numberOfNodes =
      (size * order + 1) * (size * order + 1) * (size * order + 1);
  // init new mesh
  CartesianSEMmesh<float, int, int, order> mesh(size, size, size, 200., 200.,
                                                200., order);

  // security checks
  assert(numberOfElement == mesh.getNumberOfElements());
  assert(numberOfNodes == mesh.getNumberOfNodes());

  // 1) Checking if Element To node map is valide
  vector<vector<int>> oldGlobalNodeList(numberOfElement,
                                        vector<int>(numberOfNodes, 0));
  globalNodesList(size, 2, oldGlobalNodeList);
  for (int y = 0; y < size; y++) {
    for (int z = 0; z < size; z++) {
      for (int x = 0; x < size; x++) {
        int element = x + z * size + y * size * size;
        assert(element < size * size * size);
        if (oldGlobalNodeList[element][0] !=
            mesh.globalNodeIndex(element, 0, 0, 0)) {
          cout << "Element (" << element << ") (" << x << ", " << y << ", " << z
               << ") first node: ";
          cout << mesh.globalNodeIndex(element, 0, 0, 0) << endl;
          cout << "Old implementation is: " << oldGlobalNodeList[element][0]
               << endl;
          return 1;
        };
      }
    }
  }

  // 2) Checking Coord
  // filling old data
  vector<vector<float>> oldCoodX(numberOfElement,
                                 vector<float>(numberOfNodes, 0));
  vector<vector<float>> oldCoodY(numberOfElement,
                                 vector<float>(numberOfNodes, 0));
  vector<vector<float>> oldCoodZ(numberOfElement,
                                 vector<float>(numberOfNodes, 0));
  nodesCoordinates(oldCoodX, oldCoodZ, oldCoodY, size, order, elemSize);

  // Checking X
  for (int y = 0; y < size; y++) {
    for (int z = 0; z < size; z++) {
      for (int x = 0; x < size; x++) {
        int elemIndex = x + z * size + y * size * size;
        for (int ny = 0; ny < order + 1; ny++) {
          for (int nz = 0; nz < order + 1; nz++) {
            for (int nx = 0; nx < order + 1; nx++) {
              int localNodeIndex =
                  nx + nz * (order + 1) + (order + 1) * (order + 1) * ny;
              int globalNodeIndex = mesh.globalNodeIndex(elemIndex, nx, ny, nz);

              float oldCoordxN = oldCoodX[elemIndex][localNodeIndex];
              float newCoordxN = mesh.nodeCoordX(globalNodeIndex);

              if (oldCoordxN != newCoordxN) {
                cout << "Node global id (" << globalNodeIndex << ")" << endl;
                cout << "Element (" << elemIndex << ") (" << nx << ", " << ny
                     << ", " << nz << ") node X coordinate: ";
                cout << newCoordxN << endl;
                cout << "Old implementation is: " << oldCoordxN << endl;
                return 1;
              }
            }
          }
        }
      }
    }
  }

  // Checking Y
  for (int y = 0; y < size; y++) {
    for (int z = 0; z < size; z++) {
      for (int x = 0; x < size; x++) {
        int elemIndex = x + z * size + y * size * size;
        for (int ny = 0; ny < order + 1; ny++) {
          for (int nz = 0; nz < order + 1; nz++) {
            for (int nx = 0; nx < order + 1; nx++) {
              int localNodeIndex =
                  nx + nz * (order + 1) + (order + 1) * (order + 1) * ny;
              int globalNodeIndex = mesh.globalNodeIndex(elemIndex, nx, ny, nz);

              float oldCoordyN = oldCoodY[elemIndex][localNodeIndex];
              float newCoordyN = mesh.nodeCoordY(globalNodeIndex);

              if (oldCoordyN != newCoordyN) {
                cout << "Node global id (" << globalNodeIndex << ")" << endl;
                cout << "Element (" << elemIndex << ") (" << nx << ", " << ny
                     << ", " << nz << ") node Y coordinate: ";
                cout << newCoordyN << endl;
                cout << "Old implementation is: " << oldCoordyN << endl;
                cout << "WHILE CHECKING Y" << endl;
                return 1;
              }
            }
          }
        }
      }
    }
  }

  // Checking Z
  for (int y = 0; y < size; y++) {
    for (int z = 0; z < size; z++) {
      for (int x = 0; x < size; x++) {
        int elemIndex = x + z * size + y * size * size;
        for (int ny = 0; ny < order + 1; ny++) {
          for (int nz = 0; nz < order + 1; nz++) {
            for (int nx = 0; nx < order + 1; nx++) {
              int localNodeIndex =
                  nx + nz * (order + 1) + (order + 1) * (order + 1) * ny;
              int globalNodeIndex = mesh.globalNodeIndex(elemIndex, nx, ny, nz);

              float oldCoordzN = oldCoodZ[elemIndex][localNodeIndex];
              float newCoordzN = mesh.nodeCoordZ(globalNodeIndex);

              if (oldCoordzN != newCoordzN) {
                cout << "Node global id (" << globalNodeIndex << ")" << endl;
                cout << "Element (" << elemIndex << ") (" << nx << ", " << ny
                     << ", " << nz << ") node Z coordinate: ";
                cout << newCoordzN << endl;
                cout << "Old implementation is: " << oldCoordzN << endl;
                cout << "WHILE CHECKING Z" << endl;
                return 1;
              }
            }
          }
        }
      }
    }
  }

  return 0;
}
