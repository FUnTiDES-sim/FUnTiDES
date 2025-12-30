#ifndef SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_
#define SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_

namespace model
{
template <typename Coord, typename Index>
struct CartesianParams
{
  int order;
  Index ex, ey, ez;
  Coord lx, ly, lz;
  bool isModelOnNodes;
  bool isElastic;

  Coord origin_x{0}, origin_y{0}, origin_z{0};

  CartesianParams() = default;

  CartesianParams(int order_, Index ex_, Index ey_, Index ez_, Coord lx_,
                  Coord ly_, Coord lz_, bool isModelOnNodes_, bool isElastic_)
      : order(order_),
        ex(ex_),
        ey(ey_),
        ez(ez_),
        lx(lx_),
        ly(ly_),
        lz(lz_),
        isModelOnNodes(isModelOnNodes_),
        isElastic(isElastic_)

  {
  }
};
}  // namespace model
#endif  // SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_
