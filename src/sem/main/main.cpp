//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include "SEMproxy.hpp"
#ifdef USE_EZV
#include <ezv/ezv.h>
#include <ezv/ezv_event.h>
#include <vector>
#endif // USE_EZV

#ifdef USE_CALIPER
#include "caliperUtils.hpp"
#include <caliper/cali-manager.h>
#endif // USE_CALIPER

#ifdef USE_EZV
const unsigned int SCR_WIDTH = 1024;
const unsigned int SCR_HEIGHT = 768;
static unsigned vi = 0;
static unsigned ti = 0;
static unsigned nbt = 0;
static unsigned nbv = 0;

static ezv_ctx_t ctx[2] = {NULL, NULL};
static unsigned nb_ctx = 1;
static int hud = -1;

static void do_pick(void) {
  int p = ezv_perform_1D_picking(ctx, nb_ctx);
  if (p == -1)
    ezv_hud_off(ctx[0], hud);
  else {
    ezv_hud_on(ctx[0], hud);
    ezv_hud_set(ctx[0], hud, "Cell: %d", p);
  }
}

static void process_events(void) {
  SDL_Event event;
  int r = ezv_get_event(&event, 1);
  if (r > 0) {
    int pick;
    ezv_process_event(ctx, nb_ctx, &event, NULL, &pick);
    if (pick)
      do_pick();
  }
}

struct Vec3 {
  float x, y, z;
  bool operator==(const Vec3 &other) const {
    return x == other.x && y == other.y && z == other.z;
  }
};

// Hash function for Vec3 to store unique vertices
namespace std {
template <> struct hash<Vec3> {
  size_t operator()(const Vec3 &v) const {
    return hash<float>()(v.x) ^ hash<float>()(v.y) ^ hash<float>()(v.z);
  }
};
} // namespace std

static int add_vertice(mesh3d_obj_t *mesh, float x, float y, float z) {
  mesh->vertices[vi++] = x;
  mesh->vertices[vi++] = y;
  mesh->vertices[vi++] = z;
  return nbv++;
}

static int add_triangle(mesh3d_obj_t *mesh, unsigned v1, unsigned v2,
                        unsigned v3) {
  mesh->triangles[ti++] = v1;
  mesh->triangles[ti++] = v2;
  mesh->triangles[ti++] = v3;
  return nbt++;
}

void convertSEMToMesh3D(const SEMmesh &mesh, mesh3d_obj_t *mesh3d,
                        vector<Vec3> &vertices, vector<int> &triangles) {
  unordered_map<Vec3, int> vertexIndex; // Map to store unique vertices
  int nx = mesh.getNx(), ny = mesh.getNy(), nz = mesh.getNz();
  float hx = mesh.getDx(), hy = mesh.getDy(), hz = mesh.getDz();

  mesh3d->nb_vertices = nx * ny * nz * 8;
  mesh3d->vertices = (float *)malloc(mesh3d->nb_vertices * 3 * sizeof(float));
  mesh3d->nb_triangles = nx * ny * nz * 12;
  mesh3d->triangles =
      (unsigned int *)malloc(mesh3d->nb_triangles * 3 * sizeof(unsigned));
  mesh3d->triangle_info =
      (unsigned int *)calloc(mesh3d->nb_triangles, sizeof(unsigned));

  // Iterate over all elements
  for (int i = 0; i < nx - 1; i++) {
    for (int j = 0; j < ny - 1; j++) {
      for (int k = 0; k < nz - 1; k++) {
        // Define the 8 vertices of the cube
        Vec3 cubeVertices[8] = {{i * hx, j * hy, k * hz},
                                {(i + 1) * hx, j * hy, k * hz},
                                {(i + 1) * hx, (j + 1) * hy, k * hz},
                                {i * hx, (j + 1) * hy, k * hz},
                                {i * hx, j * hy, (k + 1) * hz},
                                {(i + 1) * hx, j * hy, (k + 1) * hz},
                                {(i + 1) * hx, (j + 1) * hy, (k + 1) * hz},
                                {i * hx, (j + 1) * hy, (k + 1) * hz}};

        int vertexIndices[8];

        // Store unique vertices
        for (int v = 0; v < 8; v++) {
          if (vertexIndex.find(cubeVertices[v]) == vertexIndex.end()) {
            vertexIndex[cubeVertices[v]] = vertices.size();
            vertices.push_back(cubeVertices[v]);
            add_vertice(mesh3d, cubeVertices[v].x, cubeVertices[v].y,
                        cubeVertices[v].z);
          }
          vertexIndices[v] = vertexIndex[cubeVertices[v]];
        }

        // Define the 12 triangles (each face split into 2 triangles)
        int faces[12][3] = {
            {vertexIndices[0], vertexIndices[1], vertexIndices[2]},
            {vertexIndices[0], vertexIndices[2], vertexIndices[3]},
            {vertexIndices[1], vertexIndices[5], vertexIndices[6]},
            {vertexIndices[1], vertexIndices[6], vertexIndices[2]},
            {vertexIndices[5], vertexIndices[4], vertexIndices[7]},
            {vertexIndices[5], vertexIndices[7], vertexIndices[6]},
            {vertexIndices[4], vertexIndices[0], vertexIndices[3]},
            {vertexIndices[4], vertexIndices[3], vertexIndices[7]},
            {vertexIndices[3], vertexIndices[2], vertexIndices[6]},
            {vertexIndices[3], vertexIndices[6], vertexIndices[7]},
            {vertexIndices[4], vertexIndices[5], vertexIndices[1]},
            {vertexIndices[4], vertexIndices[1], vertexIndices[0]}};

        // Append triangles to the list
        for (int t = 0; t < 12; t++) {
          triangles.push_back(faces[t][0]);
          triangles.push_back(faces[t][1]);
          triangles.push_back(faces[t][2]);

          add_triangle(mesh3d, faces[t][0], faces[t][1], faces[t][2]);
        }
      }
    }
  }
}
#endif // USE_EZV

int main(int argc, char *argv[]) {

  time_point<system_clock> startInitTime = system_clock::now();

#ifdef USE_KOKKOS
  cout << "Using Kokkos" << endl;
  Kokkos::initialize(argc, argv);
  {
#endif

#ifdef USE_CALIPER
    cali::ConfigManager mgr;
    int caliperInitRet = launch_caliper_ctx(argc, argv, mgr);
    CALI_CXX_MARK_FUNCTION;
#endif // USE_CALIPER

    cout << "\n+================================= " << endl;
    cout << "| Initializing SEM Application ... " << endl;
    cout << "+================================= \n" << endl;

    SEMproxy semsim(argc, argv);

    semsim.initFiniteElem();

#ifdef USE_EZV
    cout << "Initialize ezv and its mesh." << endl;
    ezv_init();
    static mesh3d_obj_t mesh;

    // mesh3d_obj_build_torus_volume(&mesh, 32, 16, 16);
    std::vector<Vec3> vertices;
    std::vector<int> triangles;
    convertSEMToMesh3D(semsim.myMesh, &mesh, vertices, triangles);

    // Create SDL windows and initialize OpenGL context
    ctx[0] = ezv_ctx_create(EZV_CTX_TYPE_MESH3D, "Mesh", SDL_WINDOWPOS_CENTERED,
                            SDL_WINDOWPOS_UNDEFINED, SCR_WIDTH, SCR_HEIGHT,
                            EZV_ENABLE_PICKING | EZV_ENABLE_HUD |
                                EZV_ENABLE_CLIPPING);
    hud = ezv_hud_alloc(ctx[0]);
    ezv_hud_on(ctx[0], hud);
    // Attach mesh
    ezv_mesh3d_set_mesh(ctx[0], &mesh);
    ezv_use_data_colors_predefined(ctx[0], EZV_PALETTE_RAINBOW);

    while (1) {
      process_events();
      ezv_render(ctx, nb_ctx);
    }
    cout << "End of EZV init." << endl;
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
