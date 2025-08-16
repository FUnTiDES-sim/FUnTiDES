#ifndef SRC_MODEL_CARTESIANSEMMESH_INCLUDE_GLLPOINTS_HPP_
#define SRC_MODEL_CARTESIANSEMMESH_INCLUDE_GLLPOINTS_HPP_

#include <array>
#include <cstddef>

constexpr int MAX_GLL_ORDER = 5;

template<int Order>
struct GLLPoints {
    static_assert(Order >= 1 && Order <= MAX_GLL_ORDER,
                  "Order must be between 1 and 5");
    static constexpr int num_points = Order + 1;
    static constexpr float sqrt5 = 2.2360679774997897;

    /// Return the i-th GLL point for given Order
    PROXY_HOST_DEVICE static constexpr float get(int i) {
        if constexpr (Order == 1) {
            switch (i) {
                case 0: return -1.0;
                case 1: return 1.0;
            }
        } else if constexpr (Order == 2) {
            switch (i) {
                case 0: return -1.0;
                case 1: return 0.0;
                case 2: return 1.0;
            }
        } else if constexpr (Order == 3) {
            switch (i) {
                case 0: return -1.0;
                case 1: return -1.0 / sqrt5;
                case 2: return 1.0 / sqrt5;
                case 3: return 1.0;
            }
        } else if constexpr (Order == 4) {
            switch (i) {
                case 0: return -1.0;
                case 1: return -0.654653670707977;
                case 2: return 0.0;
                case 3: return 0.654653670707977;
                case 4: return 1.0;
            }
        } else if constexpr (Order == 5) {
            switch (i) {
                case 0: return -1.0;
                case 1: return -0.765055323929465;
                case 2: return -0.285231516480645;
                case 3: return 0.285231516480645;
                case 4: return 0.765055323929465;
                case 5: return 1.0;
            }
        }
        return 0.0;  // fallback for invalid i
    }
};

#endif  // SRC_MODEL_CARTESIANSEMMESH_INCLUDE_GLLPOINTS_HPP_
