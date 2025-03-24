//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include "SEMproxy.hpp"
#ifdef USE_CALIPER
#include <caliper/cali-manager.h>
#include <caliper/cali.h>
#include <string>
#endif // USE_CALIPER
#ifdef USE_EZV
#include <ezv/ezv.h>
#include <ezv/ezv_event.h>
#endif // USE_EZV

#ifdef USE_EZV
const unsigned int SCR_WIDTH = 1024;
const unsigned int SCR_HEIGHT = 768;

void convertSEMmeshToEZMmesh(SEMmesh semMesh, mesh3d_obj_t *mesh) { return; }
#endif // USE_EZV

int main(int argc, char *argv[]) {

  time_point<system_clock> startInitTime = system_clock::now();

#ifdef USE_KOKKOS
  cout << "Using Kokkos" << endl;
  Kokkos::initialize(argc, argv);
  {
#endif

#ifdef USE_CALIPER
    // Defines CALIPER configuration with cmdline arguments
    // Sets -P option to choose caliper outputs
    cali::ConfigManager mgr;

    std::string cali_configuration;
    if (cmdOptionExists(argv, argv + argc, "-P")) {
      cali_configuration = getCmdOption(argv, argc + argv, "-P");
    } else {
      cali_configuration = "runtime-report";
    }

    mgr.add(cali_configuration.c_str());
    if (mgr.error()) {
      std::cerr << "Config error: " << mgr.error_msg() << std::endl;
    } else {
      std::cout << "Starting Caliper with option set at: " << cali_configuration
                << std::endl;
    }
    mgr.start();
    CALI_CXX_MARK_FUNCTION;
#endif // USE_CALIPER

    cout << "\n+================================= " << endl;
    cout << "| Initializing SEM Application ... " << endl;
    cout << "+================================= \n" << endl;

    SEMproxy semsim(argc, argv);

    semsim.initFiniteElem();

#ifdef USE_EZV
    ezv_init();
    static mesh3d_obj_t mesh;
    static ezv_ctx_t ctx[2] = {NULL, NULL};
    static unsigned nb_ctx = 1;
    static int hud = -1;

    mesh3d_obj_build_torus_volume(&mesh, 32, 16, 16);
    convertSEMmeshToEZMmesh(semsim.myMesh, &mesh);

      // Create SDL windows and initialize OpenGL context
    ctx[0] = ezv_ctx_create (EZV_CTX_TYPE_MESH3D, "Mesh", SDL_WINDOWPOS_CENTERED,
			     SDL_WINDOWPOS_UNDEFINED, SCR_WIDTH, SCR_HEIGHT,
			     EZV_ENABLE_PICKING | EZV_ENABLE_HUD |
			     EZV_ENABLE_CLIPPING);
    hud    = ezv_hud_alloc (ctx[0]);
    ezv_hud_on (ctx[0], hud);
    // Attach mesh
    ezv_mesh3d_set_mesh (ctx[0], &mesh);
    ezv_use_data_colors_predefined (ctx[0], EZV_PALETTE_RAINBOW);
    ezv_render (ctx, nb_ctx);
#endif // USE_EZV

    cout << "\n+================================= " << endl;
    cout << "| Running SEM Application ...      " << endl;
    cout << "+================================= \n" << endl;

    // start timer
    time_point<system_clock> startRunTime = system_clock::now();
    semsim.run();

    cout << "\n+================================= " << endl;
    cout << "| SEM Application Finished.       " << endl;
    cout << "+================================= \n" << endl;

    // print timing information
    cout << "Elapsed Initial Time : "
         << (startRunTime - startInitTime).count() / 1E9 << " seconds." << endl;
    cout << "Elapsed Compute Time : "
         << (system_clock::now() - startRunTime).count() / 1E9 << " seconds."
         << endl;

#ifdef USE_CALIPER
    mgr.flush();
#endif // USE_CALIPER

#ifdef USE_KOKKOS
  }
  Kokkos::finalize();
#endif

  cout << "Elapsed TotalExe Time : "
       << (system_clock::now() - startInitTime).count() / 1E9 << " seconds.\n"
       << endl;
  return (0);
}
