#ifndef SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_
#define SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_

template <typename Coord, typename Index>
struct CartesianUnstructParams
{
    int order;

    Index ex, ey, ez;
    Coord lx, ly, lz;
};

#endif  // SRC_MODEL_CARTESIANUNSTRUCTMESH_INCLUDE_CARTESIAN_UNSTRUCT_PARAMS_H_
