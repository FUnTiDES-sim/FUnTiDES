#pragma once

template <typename Coord, typename Index>
struct CartesianParams :
    public BaseMesh<Coord, Index>::DataStruct
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
