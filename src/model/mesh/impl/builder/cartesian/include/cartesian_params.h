#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARAMS_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARAMS_H_
namespace model {
template <typename Coord, typename Index>
struct CartesianParams {
  int order;
  Index ex, ey, ez;
  Coord lx, ly, lz;
  bool isModelOnNodes;
  bool isElastic;
  std::string model_file;
  bool isAcoustoElastic{false};
  Coord acoustoElasticBoundaryZ{static_cast<Coord>(0)};
  Coord DgSemBoundaryZ{static_cast<Coord>(0)};

  // Global domain info (for MPI decomposition)
  Coord global_lx{0}, global_ly{0}, global_lz{0};
  Coord global_origin_x{0}, global_origin_y{0}, global_origin_z{0};

  // Local subdomain origin
  Coord origin_x{0}, origin_y{0}, origin_z{0};

  CartesianParams() = default;

  CartesianParams(int order_, Index ex_, Index ey_, Index ez_, Coord lx_, Coord ly_, Coord lz_, bool isModelOnNodes_,
                  bool isElastic_)
      : order(order_),
        ex(ex_),
        ey(ey_),
        ez(ez_),
        lx(lx_),
        ly(ly_),
        lz(lz_),
        isModelOnNodes(isModelOnNodes_),
        isElastic(isElastic_),
        model_file{""} {}
};
}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_PARAMS_H_
