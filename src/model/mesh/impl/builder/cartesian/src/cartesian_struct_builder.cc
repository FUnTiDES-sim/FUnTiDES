/**
 * @file cartesian_struct_builder.cc
 * @brief Explicit instantiations of model::CartesianStructBuilder for float/int and orders 1 to 9.
 */
#include "cartesian_struct_builder.h"

template class model::CartesianStructBuilder<float, int, 1>;
template class model::CartesianStructBuilder<float, int, 2>;
template class model::CartesianStructBuilder<float, int, 3>;
template class model::CartesianStructBuilder<float, int, 4>;
template class model::CartesianStructBuilder<float, int, 5>;
template class model::CartesianStructBuilder<float, int, 6>;
template class model::CartesianStructBuilder<float, int, 7>;
template class model::CartesianStructBuilder<float, int, 8>;
template class model::CartesianStructBuilder<float, int, 9>;
