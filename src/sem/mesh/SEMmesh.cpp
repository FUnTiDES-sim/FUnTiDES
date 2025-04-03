#include "SEMmesh.hpp"

SEMmesh::SEMmesh( const int & ex_in, const int & ey_in, const int & ez_in,
                  const float & lx_in, const float & ly_in, const float & lz_in,
                  const int & order_in )
{
  orderx=order_in;
  ordery=order_in;
  orderz=order_in;
  order=order_in;
  ex=ex_in;
  ey=ey_in;
  ez=ez_in;
  lx=lx_in;
  ly=ly_in;
  lz=lz_in;
  nx=ex_in*orderx+1;
  ny=((ey==0)?1:ey_in*ordery+1);
  nz=ez_in*orderz+1;
  hx=lx/(1.*ex);
  hy=((ey==0)?1:ly/(1.*ey));
  hz=lz/(1.*ez);
  cout <<"SEMmesh initiliazed\n";
  cout<<"ex, ey, ez="<<ex<<", "<<ey<<", "<<ez<<endl;
  cout<<"nx, ny, nz="<<nx<<", "<<ny<<", "<<nz<<endl;
}

int SEMmesh::getNumberOfNodes() const
{
  int numberOfNodes=(ex*orderx+1)*(ey*ordery+1)*(ez*orderz+1);
  printf( "Number of nodes: %d\n", numberOfNodes );
  return numberOfNodes;
}

int SEMmesh::getNumberOfElements() const
{
  int numberOfElements=((ey==0)?ex*ez:ex*ey*ez);
  printf( "Number of elements: %d\n", numberOfElements );
  return numberOfElements;
}

// get number of points per element
int SEMmesh::getNumberOfPointsPerElement() const
{
  return(((ny==1)?(orderx+1)*(orderz+1):(orderx+1)*(ordery+1)*(orderz+1)));
}

//get number of interior elements
int SEMmesh::getNumberOfInteriorElements() const
{return ((ny==1)?(ex-2)*(ez-2):(ex-2)*(ey-2)*(ez-2));}

//get number of interior Nodes
int SEMmesh::getNumberOfInteriorNodes() const
{return ((ny==1)?(nx-2)*(nz-2):(nx-2)*(ny-2)*(nz-2));}

//get nx
int SEMmesh::getNx() const
{return nx;}

//get ny
int SEMmesh::getNy() const
{return ny;}

//get nz
int SEMmesh::getNz() const
{return nz;}

//get dx
int SEMmesh::getDx() const
{return hx;}

//get dy
int SEMmesh::getDy() const
{return hy;}

//get dz
int SEMmesh::getDz() const
{return hz;}

//get coord in one direction
std::vector< float > SEMmesh::getCoordInOneDirection( const int & order, const int & nCoord, const int & h, const int & nElement ) const
{
  std::vector< float > coord( nCoord );
  std::vector< float > xi( order+1 );
  switch( order )
  {
    case 1:
      xi[0]=-1.;
      xi[1]=1.;
      break;
    case 2:
      xi[0]=-1.;
      xi[1]=0.;
      xi[2]=1.;
      break;
    case 3:
      static constexpr double sqrt5 = 2.2360679774997897;
      xi[0] = -1.0;
      xi[1] = -1./sqrt5;
      xi[2] = 1./sqrt5;
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
      xi[1] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
      xi[2] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
      xi[3] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
      xi[4] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
      xi[5] = 1.0;
      break;
    default:
      break;
  }
  int i=nElement;
  float x0=i*h;
  float x1=(i+1)*h;
  float b=(x1+x0)/2.;
  float a=b-x0;
  for( int j=0; j<order+1; j++ )
  {
    coord[j]=a*xi[j]+b;
  }
  return coord;
}

/*
// Initialize nodal coordinates.
void SEMmesh::nodesCoordinates( arrayReal & nodeCoordsX,
                                arrayReal & nodeCoordsZ,
                                arrayReal & nodeCoordsY) const
{
  std::vector< float > coordX( order+1 );
  std::vector< float > coordY( order+1 );
  std::vector< float > coordZ( order+1 );

  for(int n=0;n<ey;n++)
  {
     coordY=getCoordInOneDirection( order, order+1, hy, n );
     for(int m=0;m<ez;m++)
     {
        coordZ=getCoordInOneDirection( order, order+1, hz, m );
        for(int l=0;l<ex;l++)
        {
           coordX=getCoordInOneDirection( order, order+1, hx, l );
           int e=l+m*ex+n*ex*ez;
           for( int k=0; k<order+1; k++ )
           {
              for( int j=0; j<order+1; j++ )
              {
                 for( int i=0; i<order+1; i++ )
                 {
                    nodeCoordsX(e, i+(order+1)*j+k*(order+1)*(order+1))=coordX[i];
                    nodeCoordsZ(e, i+(order+1)*j+k*(order+1)*(order+1))=coordZ[j];
                    nodeCoordsY(e, i+(order+1)*j+k*(order+1)*(order+1))=coordY[k];
                 }
              }
           }
        }
     }
  }
}
*/
// Initialize nodal coordinates.
void SEMmesh::nodesCoordinates(arrayReal & nodeCoordsX,
                               arrayReal & nodeCoordsZ,
                               arrayReal & nodeCoordsY) const
{

  for(int n=0;n<ey;n++)
  {
    for(int m=0;m<ez; m++)
    {
      for(int l=0;l<ex;l++)
      {
        int e=l+m*ex+n*ex*ez;
        for( int k=0; k<2; k++ )
        {
          for( int j=0; j<2; j++ )
          {
            for( int i=0; i<2; i++ )
            {
              nodeCoordsX(e, i+2*j+4*k)=(l+i)*hx;
              nodeCoordsZ(e, i+2*j+4*k)=(m+j)*hz;
              nodeCoordsY(e, i+2*j+4*k)=(n+k)*hy;
            }
          }
        }
      }
    }
  }
}

//  list of global nodes ( vertices) for each element
void SEMmesh::globalNodesList( const int & numberOfElements, arrayInt & nodesList ) const
{
  for( int j=0; j<((ey==0)?1:ey); j++ )
  {
    for( int k=0; k<ez; k++ )
    {
      for( int i=0; i<ex; i++ )
      {
        int n0=i+k*ex+j*ex*ez;
        int offset=i*order+k*order*nx+j*order*nx*nz;
        for( int m=0; m<((ey==0)?1:order+1); m++ )
        {
          for( int n=0; n<order+1; n++ )
          {
            for( int l=0; l<order+1; l++ )
            {
              int dofLocal=l+n*(order+1)+m*(order+1)*(order+1);
              int dofGlobal=offset+l+n*nx+m*nx*nz;
              nodesList( n0, dofLocal )=dofGlobal;
              //if(n0==7)printf("offset=%d dofLocal=%d dofGlobal=%d\n",offset,dofLocal,dofGlobal);
            }
          }
        }
      }
    }
  }
}


// compute global node to grid  indexes
int SEMmesh::Itoijk( const int & I, int & i, int & j, int & k ) const
{
  i=I%nx;
  j=int((I%(nx*nz))/nx);
  k=I/(nx*nz);
  return 0;
}

// compute global node to grid  indexes
int SEMmesh::Itoij( const int & I, int & i, int & j ) const
{
  i=I%nx;
  j=int((I-i)/nx);
  return 0;
}

// project vector node to grid
std::vector< std::vector< float > > SEMmesh::projectToGrid( const int numberOfNodes, const std::vector< float > inputVector ) const
{
  std::vector< vector< float > > grid( nx, std::vector< float >( nz ));
  int i, j;
  //cout<< " In projecToGrid numberOfNodes="<<numberOfNodes<<endl;
  for( int node=0; node<numberOfNodes; node++ )
  {
    Itoij( node, i, j );
    grid[i][j]=inputVector[node];
  }
  return grid;
}

// compute element e where (x,y,z) belongs to
int SEMmesh::getElementNumberFromPoints( const float & x, const float & y, const float & z ) const
{
  int eX=0, eY=0, eZ=0;
  for( int i=0; i<ex; i++ )
  {
    if( x>= i*hx && x<(i+1)*hx )
      eX=i;
  }
  for( int i=0; i<ey; i++ )
  {
    if( y>= i*hy && y<(i+1)*hy )
      eY=i;
  }
  for( int i=0; i<ez; i++ )
  {
    if( z>= i*hz && z<(i+1)*hz )
      eZ=i;
  }
  return eX+eZ*ex+eY*ex*ez;
}

// set model
void SEMmesh::getModel( const int & numberOfElements, vectorReal & model ) const
{

  for( int j=0; j<((ey==0)?1:ey); j++ )
  {
    for( int k=0; k<ez; k++ )
    {
      for( int i=0; i<ex; i++ )
      {
        int e=i+k*ex+j*ex*ez;
        model[e]=1500;
      }
    }
    /*
       for( int k=ez/2; k<ez; k++ )
       {
       for( int i=0; i<ex; i++ )
       {
         int e=i+k*ex+j*ex*ez;
         model[e]=3500;
       }
       }
     */
  }
}
//  get list of global interior nodes
int SEMmesh::getListOfInteriorNodes( const int & numberOfInteriorNodes,
                                     vectorInt & listOfInteriorNodes ) const
{
  int m=0;
  if( ny==1 )
  {
    for( int j=1; j<nz-1; j++ )
    {
      for( int i=1; i<nx-1; i++ )
      {
        listOfInteriorNodes[m]=i+j*nx;
        m++;
      }
    }
  }
  else
  {
    for( int k=1; k<ny-1; k++ )
    {
      for( int j=1; j<nz-1; j++ )
      {
        for( int i=1; i<nx-1; i++ )
        {
          listOfInteriorNodes[m]=i+j*nx+k*nx*nz;
          m++;
        }
      }
    }
  }
  return 0;
}

// get list of interior Elements
void SEMmesh::getListOfInteriorElements( vectorInt & listOfInteriorElements ) const
{
  int m=0;
  if( ey==0 )
  {
    for( int j=1; j<ez-1; j++ )
    {
      for( int i=1; i<ex-1; i++ )
      {
        listOfInteriorElements[m]=i+j*ex;
        m++;
      }
    }
  }
  else
  {
    for( int k=1; k<ey-1; k++ )
    {
      for( int j=1; j<ez-1; j++ )
      {
        for( int i=1; i<ex-1; i++ )
        {
          listOfInteriorElements[m]=i+j*ex+k*ex*ez;
          m++;
        }
      }
    }
  }
}


// save snapshot
void SEMmesh::saveSnapShot( const int indexTimeStep, const int i1, arrayReal const & u ) const
{

  int nx=getNx();
  int ny=getNy();
  int nz=getNz();
  float dx=getDx();
  float dy=getDy();
  float dz=getDz();
  int numberOfNodes=nx*nz;
  std::vector< float > inputVector( numberOfNodes );
  int offset=(ny==1?offset=0:offset=nx*nz*(ny/2-1));
  //printf(" nx, ny/2-1, nz %d %d %d offset=%d\n",nx,ny/2-1,nz, offset);
  for( int i = offset; i< offset+numberOfNodes; i++ )
  {
    inputVector[i-offset]=u( i, i1 );
  }
  std::vector< std::vector< float > > grid=projectToGrid( numberOfNodes, inputVector );
  fstream snapFile;
  string snapNumber = "snapshot"+to_string( indexTimeStep );
  snapFile.open( snapNumber, ios::out| ios::trunc );
  //std::cout<<"nx="<<nx<<" ny="<<ny<<std::endl;
  for( int i=0; i<nx; i++ )
  {
    for( int j=0; j<nz; j++ )
    {
      snapFile<<i*dx<<" "<<j*dz<<" " <<grid[i][j]<<endl;
    }
  }
  snapFile.close();
}
