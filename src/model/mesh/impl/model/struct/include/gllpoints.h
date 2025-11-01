#ifndef SRC_MODEL_MODELAPI_INCLUDE_GLLPOINTS_H_
#define SRC_MODEL_MODELAPI_INCLUDE_GLLPOINTS_H_

#include <array>
#include <cstddef>

constexpr int MAX_GLL_ORDER = 5;

struct GLLPoints
{
  static constexpr float sqrt5 = 2.2360679774997897;

  /// Return the number of points for given order
  PROXY_HOST_DEVICE static constexpr int num_points(int order)
  {
    return order + 1;
  }

  /// Return the i-th GLL point for given order
  PROXY_HOST_DEVICE static float get(int order, int i)
  {
    // Validate order
    if (order < 1 || order > MAX_GLL_ORDER)
    {
      return 0.0f;  // Invalid order
    }

    // Validate point index
    if (i < 0 || i >= (order + 1))
    {
      return 0.0f;  // Invalid index
    }

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
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -1.0f / sqrt5;
          case 2:
            return 1.0f / sqrt5;
          case 3:
            return 1.0f;
        }
        break;

      case 4:
        switch (i)
        {
          case 0:
            return -1.0f;
          case 1:
            return -0.654653670707977f;
          case 2:
            return 0.0f;
          case 3:
            return 0.654653670707977f;
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
            return -0.765055323929465f;
          case 2:
            return -0.285231516480645f;
          case 3:
            return 0.285231516480645f;
          case 4:
            return 0.765055323929465f;
          case 5:
            return 1.0f;
        }
        break;
    }

    return 0.0f;  // fallback for invalid cases
  }
};

#endif  // SRC_MODEL_MODELAPI_INCLUDE_GLLPOINTS_H_
