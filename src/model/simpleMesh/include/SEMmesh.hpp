#ifndef SEM_MESH_
#define SEM_MESH_

#include <SEMmacros.hpp>
#include <dataType.hpp>
using namespace std;

class SEMmesh {
private:
  int ex, ey, ez;
  int nx, ny, nz;
  float lx, ly, lz;
  float hx, hy, hz;
  int orderx, ordery, orderz, order;
  int nbFaces;
  int spongeSize;
  bool surfaceSponge;
  bool surfaceDamping;

public:
  PROXY_HOST_DEVICE SEMmesh() {};
  PROXY_HOST_DEVICE ~SEMmesh(){};

  SEMmesh(const int &ex_in, const int &ey_in, const int &ez_in,
          const float &lx_in, const float &ly_in, const float &lz_in,
          const int &order_in);

  SEMmesh(const int &ex_in, const int &ey_in, const int &ez_in,
          const float &lx_in, const float &ly_in, const float &lz_in,
          const int &order_in, const int spongeSize, const bool surfaceSponge,
          const bool surfaceDamping)
      : SEMmesh(ex_in, ey_in, ez_in, lx_in, ly_in, lz_in,
                order_in) // Delegation happens here
  {
    this->spongeSize = spongeSize;
    this->surfaceSponge = surfaceSponge;
    this->surfaceDamping = surfaceDamping;
  }

  // Returns number of Nodes in the mesh
  int getNumberOfNodes() const;

  //  Returns the number of elements of the mesh
  int getNumberOfElements() const;

  // get number of points per element
  int getNumberOfPointsPerElement() const;

  // get number of interior elements
  int getNumberOfInteriorElements() const;

  int getNumSurfaceElementsToExcludeFromSponge() const;

  // get number of interior nodes
  int getNumberOfInteriorNodes() const;
  int getNumberOfInteriorNodes(int spongeSize) const;

  // get number of sponge elements
  int getNumberOfSpongeElements() const;

  // get number of sponge nodes
  int getNumberOfSpongeNodes() const;

  // get number of damping nodes
  // TODO
  int getNumberOfDampingNodes() const;

  // get number of damping Elements
  int getNumberOfDampingElements() const;

  // get number of bundary nodes
  int getNumberOfBundaryNodes() const;

  // get nx
  int getNx() const;

  // get ny
  int getNy() const;

  // get nz
  int getNz() const;

  // get Dx
  int getDx() const;

  // get Dy
  int getDy() const;

  // get Dz
  int getDz() const;

  // getter for ex
  int getEx() const;

  // getter for ey
  int getEy() const;

  // getter for ez
  int getEz() const;

  // get coord in one direction
  vector<float> getCoordInOneDirection(const int &order, const int &nCoord,
                                       const int &h, const int &nElement) const;

  // Initialize nodal coordinates.
  void nodesCoordinates(arrayReal &nodeCoordsX, arrayReal &nodeCoordsZ,
                        arrayReal &nodeCoordsY) const;

  //  list of global nodes ( vertices)
  void globalNodesList(const int &numberOfElements, arrayInt &nodesList) const;

  // compute element e where (x,y,z) belongs to
  int getElementNumberFromPoints(const float &x, const float &y,
                                 const float &z) const;

  // get list of interior Elements
  void getListOfInteriorElements(vectorInt &listOfInteriorElements) const;

  //  get list of global interior nodes
  int getListOfInteriorNodes(const int &numberOfInteriorNodes,
                             vectorInt &listOfInteriorNodes) const;

  // set model
  void getModel(const int &numberOfNodes, vectorReal &model) const;

  // compute global to local node index 3D version
  int Itoijk(const int &I, int &i, int &j, int &k) const;

  // compute global to local node index 2D version
  int Itoij(const int &I, int &i, int &j) const;

  // project vector node to grid
  vector<vector<float>> projectToGrid(const int numberOfNodes,
                                      const vector<float> inputVector) const;

  // save snapshot
  void saveSnapShot(const int indexTimeStep, const int i1,
                    arrayReal const &u) const;
};
#endif // SEM_MESH_
