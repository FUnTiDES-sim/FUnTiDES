#ifndef FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_GLLPOINTS_H_
#define FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_GLLPOINTS_H_
#include <array>
#include <cstddef>

constexpr int MAX_GLL_ORDER = 9;

struct GLLPoints
{
  /// Return the number of points for given order
  PROXY_HOST_DEVICE static constexpr int num_points(int order)
  {
    return order + 1;
  }

  /// Return the i-th GLL point for given order (reference element [-1, 1])
  PROXY_HOST_DEVICE static float get(int order, int i)
  {
    if (order < 1 || order > MAX_GLL_ORDER) return 0.0f;
    if (i < 0 || i >= (order + 1)) return 0.0f;

    switch (order)
    {
      case 1:
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return 1.0f;
        }
        break;

      case 2:
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return 0.0f;
          case 2:
            return 1.0f;
        }
        break;

      case 3:
        // interior: ±1/√5
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -0.4472135954999579f;
          case 2:
            return 0.4472135954999579f;
          case 3:
            return 1.0f;
        }
        break;

      case 4:
        // interior: ±√(3/7), 0
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -0.6546536707079772f;
          case 2:
            return 0.0f;
          case 3:
            return 0.6546536707079772f;
          case 4:
            return 1.0f;
        }
        break;

      case 5:
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -0.7650553239294647f;
          case 2:
            return -0.2852315164806451f;
          case 3:
            return 0.2852315164806451f;
          case 4:
            return 0.7650553239294647f;
          case 5:
            return 1.0f;
        }
        break;

      case 6:
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -0.8302238962785670f;
          case 2:
            return -0.4688487934707142f;
          case 3:
            return 0.0f;
          case 4:
            return 0.4688487934707142f;
          case 5:
            return 0.8302238962785670f;
          case 6:
            return 1.0f;
        }
        break;

      case 7:
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -0.8717401485096066f;
          case 2:
            return -0.5917001814331423f;
          case 3:
            return -0.2092992179024789f;
          case 4:
            return 0.2092992179024789f;
          case 5:
            return 0.5917001814331423f;
          case 6:
            return 0.8717401485096066f;
          case 7:
            return 1.0f;
        }
        break;

      case 8:
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -0.8997579954114602f;
          case 2:
            return -0.6771862795107377f;
          case 3:
            return -0.3631174638261782f;
          case 4:
            return 0.0f;
          case 5:
            return 0.3631174638261782f;
          case 6:
            return 0.6771862795107377f;
          case 7:
            return 0.8997579954114602f;
          case 8:
            return 1.0f;
        }
        break;

      case 9:
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -0.9195339081664589f;
          case 2:
            return -0.7387738651055050f;
          case 3:
            return -0.4779249498104445f;
          case 4:
            return -0.1652789576663870f;
          case 5:
            return 0.1652789576663870f;
          case 6:
            return 0.4779249498104445f;
          case 7:
            return 0.7387738651055050f;
          case 8:
            return 0.9195339081664589f;
          case 9:
            return 1.0f;
        }
        break;
    }
    return 0.0f;
  }
};

#endif  // FUNTIDES_MODEL_MESH_IMPL_MODEL_STRUCT_INCLUDE_GLLPOINTS_H_
