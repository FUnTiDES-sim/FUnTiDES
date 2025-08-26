#ifndef SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_
#define SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_

template <typename Coord, typename Index>
struct CartesianParams :
    public mesh_api::MeshBase::DataStruct
{
    int order;

    Index ex, ey, ez;
    Coord lx, ly, lz;

    CartesianParams() = default;

    CartesianParams(int order_,
                    Index ex_, Index ey_, Index ez_,
                    Coord lx_, Coord ly_, Coord lz_)
        : order(order_), ex(ex_), ey(ey_), ez(ez_),
          lx(lx_), ly(ly_), lz(lz_) {}
};

#endif  // SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_
