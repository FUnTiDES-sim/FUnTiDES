#ifdef USE_EZV
#include "SEMmesh.hpp"
#include "SEMproxy.hpp"
#include <ezv/ezv.h>
#include <ezv/ezv_event.h>
#include <vector>

const unsigned int SCR_WIDTH = 1024;
const unsigned int SCR_HEIGHT = 768;
static unsigned vi = 0;
static unsigned ti = 0;
static unsigned nbt = 0;
static unsigned nbv = 0;

static ezv_ctx_t ctx[2] = {NULL, NULL};
static unsigned nb_ctx = 1;
static int hud = -1;

/**
 * @brief Handles the pick operation.
 *
 * This function performs a 1D picking operation using `ezv_perform_1D_picking`.
 * If the picking operation returns `-1`, it turns off the HUD display.
 * Otherwise, it turns on the HUD and displays the picked cell's ID.
 */
static void do_pick(void) {
  int p = ezv_perform_1D_picking(ctx, nb_ctx);
  if (p == -1)
    ezv_hud_off(ctx[0], hud);
  else {
    ezv_hud_on(ctx[0], hud);
    ezv_hud_set(ctx[0], hud, const_cast<char *>("Cell: %d"), p);
  }
}

enum { EZV_THR_EVENT_DATA_COLORS, EZV_THR_EVENT_CPU_COLORS };

static uint32_t base_event = 0;

static void thr_push_data_colors(ezv_ctx_t ctx, void *values) {
  SDL_Event event;

  event.type = SDL_USEREVENT;
  event.user.code = base_event + EZV_THR_EVENT_DATA_COLORS;
  event.user.data1 = (void *)ctx;
  event.user.data2 = (void *)values;

  SDL_PushEvent(&event);
}

ezv_ctx_t get_ezv_ctx() { return ctx; }

/**
 * @brief Processes events from the event queue.
 *
 * This function retrieves events from the event queue and processes them using
 * `ezv_process_event`. If a pick event is detected, it calls the `do_pick`
 * function to handle the pick operation.
 */
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

/**
 * @brief Structure representing a 3D vector.
 *
 * This structure is used to represent a 3D vector with `x`, `y`, and `z`
 * coordinates.
 */
struct Vec3 {
  float x, y, z;

  /**
   * @brief Compares two Vec3 objects for equality.
   *
   * This operator compares the `x`, `y`, and `z` coordinates of two `Vec3`
   * objects and returns `true` if all the coordinates are equal, `false`
   * otherwise.
   *
   * @param other The `Vec3` object to compare with.
   *
   * @return `true` if the vectors are equal, `false` otherwise.
   */
  bool operator==(const Vec3 &other) const {
    return x == other.x && y == other.y && z == other.z;
  }

  // Opérateur de sortie pour std::cout
  friend std::ostream &operator<<(std::ostream &os, const Vec3 &v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
  }
};

/**
 * @brief Specialization of the hash function for Vec3 to store unique vertices.
 *
 * This specialization of the `std::hash` template allows the usage of `Vec3` as
 * a key in unordered containers by generating a unique hash value based on the
 * coordinates (x, y, z) of the `Vec3` object.
 *
 * @param v The Vec3 object to hash.
 *
 * @return A hash value representing the Vec3 object.
 */
// Hash function for Vec3 to store unique vertices
namespace std {
template <> struct hash<Vec3> {
  size_t operator()(const Vec3 &v) const {
    return hash<float>()(v.x) ^ hash<float>()(v.y) ^ hash<float>()(v.z);
  }
};
} // namespace std

/**
 * @brief Adds a vertex to the EZV mesh.
 *
 * This function adds a vertex with the given coordinates (x, y, z) to the
 * mesh's vertex list.
 *
 * @param mesh Pointer to the mesh3d object where the vertex will be added.
 * @param x The x-coordinate of the vertex.
 * @param y The y-coordinate of the vertex.
 * @param z The z-coordinate of the vertex.
 *
 * @return The number of vertices in the mesh after the addition.
 */
static int add_vertice(mesh3d_obj_t *mesh, float x, float y, float z) {
  mesh->vertices[vi++] = x;
  mesh->vertices[vi++] = y;
  mesh->vertices[vi++] = z;
  return nbv++;
}

/**
 * @brief Adds a triangle to the EZV mesh.
 *
 * @param mesh Pointer to the mesh3d object where the triangle will be added.
 * @param v1 Index of the first vertex of the triangle.
 * @param v2 Index of the second vertex of the triangle.
 * @param v3 Index of the third vertex of the triangle.
 *
 * @return The number of triangles in the mesh after the addition.
 */
static int add_triangle(mesh3d_obj_t *mesh, unsigned v1, unsigned v2,
                        unsigned v3) {
  mesh->triangles[ti++] = v1;
  mesh->triangles[ti++] = v2;
  mesh->triangles[ti++] = v3;
  return nbt++;
}

/**
 * @brief Convert SEMmesh surfacic mesh into an EZV-compatible mesh.
 *
 * Takes a SEMmesh and converts all its hexahedral elements into 12 unique
 * triangles, 2 per face. For each element of SEMmesh, it creates 8 unique
 * vertices, then uses a hash table to avoid duplicates (shared nodes).
 *
 * TODO: An optimization could be to avoid using a hash table and compute
 * boundary elements separately. This should allow for greater parallelism
 * opportunities.
 *
 * @param mesh SEMmesh mesh to convert
 * @param mesh3d EZV mesh output
 */
void convertSEMToMesh3D(const SEMmesh &mesh, mesh3d_obj_t *mesh3d) {
  vector<Vec3> vertices;
  unordered_map<Vec3, int> vertexIndex; // Map to store unique vertices
  int ex = mesh.getEx(), ey = mesh.getEy(), ez = mesh.getEz();
  float hx = mesh.getDx(), hy = mesh.getDy(), hz = mesh.getDz();

  // initial parameters
  mesh3d->nb_vertices = (ex + 1) * (ey + 1) * (ez + 1);
  mesh3d->vertices = (float *)malloc(mesh3d->nb_vertices * sizeof(float) * 3);
  mesh3d->nb_triangles = ex * ey * ez * 12; // TODO: Avoid duplicates here
  mesh3d->triangles =
      (unsigned int *)malloc(mesh3d->nb_triangles * 3 * sizeof(unsigned));
  mesh3d->triangle_info =
      (unsigned int *)calloc(mesh3d->nb_triangles, sizeof(unsigned));
  mesh3d->nb_cells = mesh.getNumberOfElements();
  mesh3d->cells = (unsigned *)malloc((mesh3d->nb_cells + 1) * sizeof(unsigned));

  int elem_id = 0;
  // Iterate over all elements
  for (int i = 0; i < ex; i++) {
    for (int j = 0; j < ey; j++) {
      for (int k = 0; k < ez; k++) {
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

        // Append triangles to the mesh
        // Add cell for this element
        mesh3d->cells[elem_id] = elem_id * 12;
        for (int t = 0; t < 12; t++) {
          int triangle_idx = t + 12 * elem_id;
          mesh3d->triangle_info[triangle_idx] |= (elem_id << CELLNO_SHIFT);
          add_triangle(mesh3d, faces[t][0], faces[t][1], faces[t][2]);
        }

        elem_id++;
      }
    }
  }
  std::cout << "Created " << vertices.size() << " vertices while "
            << (ex + 1) * (ey + 1) * (ez + 1) << " where expected."
            << std::endl;
  mesh3d->cells[elem_id] = ++elem_id * 12;
}

/**
 * @brief Initialize EZV, convert simulation mesh to ezv compatible mesh then
 * open SDL window displaying it.
 *
 * @param semsim SEMsim object containing the simulation mesh in SEMmesh format.
 * @param mesh New EZV mesh to display.
 */

inline void init_ezv() {
  cout << "Initialize ezv." << endl;
  ezv_init();

  // Create SDL windows and initialize OpenGL context
  ctx[0] =
      ezv_ctx_create(EZV_CTX_TYPE_MESH3D, "Mesh", SDL_WINDOWPOS_CENTERED,
                     SDL_WINDOWPOS_UNDEFINED, SCR_WIDTH, SCR_HEIGHT,
                     EZV_ENABLE_PICKING | EZV_ENABLE_HUD | EZV_ENABLE_CLIPPING);
  hud = ezv_hud_alloc(ctx[0]);
  ezv_hud_on(ctx[0], hud);

  cout << "End of EZV init." << endl;
}

inline void ezv_loop(void) {
  while (1) {
    process_events();
    ezv_render(ctx, nb_ctx);
  }
}

inline void ezv_init_mesh(SEMproxy &semsim, mesh3d_obj_t *mesh) {
  cout << "Converting SEMmesh to EZV mesh." << endl;
  // Setting mesh into ezm format
  convertSEMToMesh3D(semsim.myMesh, mesh);
  // Attach mesh
  ezv_mesh3d_set_mesh(ctx[0], mesh);
  ezv_use_data_colors_predefined(ctx[0], EZV_PALETTE_RAINBOW);
}

#endif // USE_EZV
