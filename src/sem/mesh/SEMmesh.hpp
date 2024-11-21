#ifndef SEM_MESH_
#define SEM_MESH_

#include "dataType.hpp"
#include "SEMmacros.hpp"
using namespace std;

class SEMmesh
{
private:
  int ex, ey, ez;
  int nx, ny, nz;
  float lx, ly, lz;
  float hx, hy, hz;
  int orderx, ordery, orderz, order;
  int nbFaces;

public:

  PROXY_HOST_DEVICE SEMmesh(){};
  PROXY_HOST_DEVICE ~SEMmesh(){};

  SEMmesh( const int & ex_in, const int & ey_in, const int & ez_in,
           const float & lx_in, const float & ly_in, const float & lz_in,
           const int & order_in );

  // Returns number of Nodes in the mesh
  int  getNumberOfNodes() const;

  //  Returns the number of elements of the mesh
  int  getNumberOfElements() const;

  // get number of points per element
  int getNumberOfPointsPerElement() const;

  //get number of interior elements
  int getNumberOfInteriorElements() const;

  //get number of interior elements
  int getNumberOfInteriorNodes() const;

  //get nx
  int getNx() const;

  //get ny
  int getNy() const;

  //get nz
  int getNz() const;

  //get Dx
  int getDx()const;

  //get Dy
  int getDy() const;

  //get Dz
  int getDz() const;

  //get coord in one direction
  vector< float > getCoordInOneDirection( const int & order, const int & nCoord, const int & h, const int & nElement ) const;


  // Initialize nodal coordinates.
  void nodesCoordinates( arrayReal & nodeCoordsX,arrayReal & nodeCoordsZ, arrayReal & nodeCoordsY ) const;

  //  list of global nodes ( vertices)
  void globalNodesList( const int & numberOfElements, arrayInt & nodesList ) const;

  // compute element e where (x,y,z) belongs to
  int getElementNumberFromPoints( const float & x, const float & y, const float & z ) const;

  // get list of interior Elements
  void getListOfInteriorElements( vectorInt & listOfInteriorElements ) const;

  //  get list of global interior nodes
  int getListOfInteriorNodes( const int & numberOfInteriorNodes, vectorInt & listOfInteriorNodes ) const;

  // set model
  void getModel( const int & numberOfNodes, vectorReal & model ) const;

  // sort element by color
  // red=0, green=1, blue=2, yellow=3
  int getNumberOfElementsByColor() const;
  void sortElementsByColor( int numberOfElementsByColor[], arrayInt & listOfElementsByColor ) const;

  // get number of Boundary Faces
  int getNumberOfBoundaryFaces() const;

  // get number of Boundary nodes
  int getNumberOfBoundaryNodes() const;

  // provides informations about boundary  faces:
  // element number,
  // orientation of the face
  void getBoundaryFacesInfos( arrayInt & faceInfos )const;

  //  get list of global boundary nodes
  int getListOfBoundaryNodes( const int & numberOfBoundaryNodes, vectorInt & listOfBoundaryNodes ) const;

  // provides a mapping between local node of a face and global node Face:
  void getLocalFaceNodeToGlobalFaceNode( arrayInt & localFaceNodeToGlobalFaceNode ) const;

  // compute global to local node index
  int Itoijk( const int & I, int & i, int & j, int & k ) const;

  // compute global to local node index
  int Itoij( const int & I, int & i, int & j ) const;

  // project vector node to grid
  vector< vector< float > > projectToGrid( const int numberOfNodes, const vector< float > inputVector ) const;

  // save snapshot
  void saveSnapShot( const int indexTimeStep, const int i1, arrayReal const & u ) const;
};
#endif //SEM_MESH_
