#ifndef FDTDGRIDS_HPP_
#define FDTDGRIDS_HPP_

#include "FDTDmacros.hpp"
#include "FDproxyOptions.hpp"
#include <dataType.hpp>

struct FDTDGRIDS
{
  // Grid dimensions
  int nx{0}, ny{0}, nz{0};
  float dx{0.0f}, dy{0.0f}, dz{0.0f};

  // Sponge / PML
  int ntaperx{0}, ntapery{0}, ntaperz{0};
  int ndampx{0} , ndampy{0}, ndampz{0};
  
  float hdx_2{0.0f}, hdy_2{0.0f}, hdz_2{0.0f};

  // boundary limits
  int x1{0}, x2{0}, x3{0}, x4{0}, x5{0}, x6{0};
  int y1{0}, y2{0}, y3{0}, y4{0}, y5{0}, y6{0};
  int z1{0}, z2{0}, z3{0}, z4{0}, z5{0}, z6{0};

  // grid physical properties
  vectorReal vp;  // p velolicty
  vectorReal vs;  // s velocity
  vectorReal rho; // density

  // init grid fo computation
  void initGrid(FDProxyOptions &m_opt)
  {
    // grid dimensions
    nx = m_opt.grid.nx;
    ny = m_opt.grid.ny;
    nz = m_opt.grid.nz;

    //grid sampling
    dx = m_opt.grid.dx;
    dy = m_opt.grid.dy;
    dz = m_opt.grid.dz;

    //get model info if from file
    getModelInfo( nx,ny,nz,
                  dx,dy,dz,
                  m_opt.velocity.fileModel);


    hdx_2=1./(4. * dx * dx);
    hdy_2=1./(4. * dy * dy);
    hdz_2=1./(4. * dz * dz);
 
  }
  
  //get model info if from file
  void getModelInfo(int &nx, int &ny, int &nz,
                          float &dx, float &dy, float &dz,
                          const std::string &fileModel)
  {
    if(fileModel!="")
    {
      // read model from file
      std::ifstream infile(fileModel, std::ios::in | std::ios::binary);
      if (!infile)
      {
        std::cerr << "Error opening file " << fileModel << std::endl;
        exit(1);
      }
      infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
      infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
      infile.read(reinterpret_cast<char*>(&nz), sizeof(int));
      infile.read(reinterpret_cast<char*>(&dx), sizeof(float));
      infile.read(reinterpret_cast<char*>(&dy), sizeof(float));
      infile.read(reinterpret_cast<char*>(&dz), sizeof(float));
      infile.close();
      printf("model read from file %s\n",fileModel.c_str());
    }
  }

  // init model and arrays
  void initModelArrays(FDProxyOptions &m_opt)
  { 

  int modelVolume = nx *ny * nz;

  vp   = allocateVector< vectorReal >( modelVolume, "vp" );
  
  // initialize velocity model
  float timeStep = 0.001f;
  float init_vp_value = m_opt.velocity.vmin*m_opt.velocity.vmin*timeStep*timeStep;
  printf("init_vp_value=%f\n",init_vp_value);

  initModel(init_vp_value,false);
  
}
  void initModel(float init_vp_value,bool fromFile=false)
  {

    if(fromFile==false)
    {    // initialize vp and pressure field
      #pragma omp parallel for collapse(3)
      for( int i=0; i<nx; i++ )
      {
        for( int j=0; j<ny; j++ )
        {
          for( int k=0; k<nz; k++ )
          {
            vp[nz*ny*(i) + nz*(j) + (k)]=init_vp_value;
          }
        }
      }
      #pragma omp parallel for collapse(3)
      for( int i=0; i<nx; i++ )
      {
        for( int j=0; j<ny; j++ )
        {
          for( int k=nz/2; k<nz; k++ )
          {
            vp[nz*ny*(i) + nz*(j) + (k)]=2*init_vp_value;
          }
        }
      }
    }
    /*else
    {
      // read vp from file
      std::ifstream infile("vp.bin", std::ios::in | std::ios::binary);
      if (!infile)
      {
        std::cerr << "Error opening file vp.bin" << std::endl;
        exit(1);
      }
      infile.read(reinterpret_cast<char*>(vp.data()), vp.size() * sizeof(float));
      infile.close();
    }*/
  }
  

};

#endif //FDTDGRIDS_HPP_
