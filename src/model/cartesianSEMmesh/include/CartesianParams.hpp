#pragma once

template< typename index_t, typename discretization_t >
struct CartesianParams
{
  int order;

  index_t ex, ey, ez;
  discretization_t lx, ly, lz;
};
