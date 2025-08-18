#ifndef SRC_MODEL_CARTESIANSEMMESH_INCLUDE_CARTESIANPARAMS_HPP_
#define SRC_MODEL_CARTESIANSEMMESH_INCLUDE_CARTESIANPARAMS_HPP_

template <typename Coord, typename Index>
struct CartesianParams :
    public mesh_base::BaseMesh<Coord, Index>::DataStruct {
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

#endif  // SRC_MODEL_CARTESIANSEMMESH_INCLUDE_CARTESIANPARAMS_HPP_
