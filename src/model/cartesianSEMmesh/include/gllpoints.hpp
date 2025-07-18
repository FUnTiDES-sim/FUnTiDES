#ifndef GLLPOINTS_H_
#define GLLPOINTS_H_

#pragma once
#include <array>
#include <cstddef>

constexpr int MAX_GLL_ORDER = 5;

template<int ORDER>
struct GLLPoints {
    static_assert(ORDER >= 1 && ORDER <= MAX_GLL_ORDER, "Order must be between 1 and 5");
    static constexpr int num_points = ORDER + 1;

    /// Return the i-th GLL point for given Order
    PROXY_HOST_DEVICE static constexpr double get(int i) {
        if constexpr (ORDER == 1) {
            switch(i) {
                case 0: return -1.0;
                case 1: return 1.0;
            }
        } else if constexpr (ORDER == 2) {
            switch(i) {
                case 0: return -1.0;
                case 1: return 0.0;
                case 2: return 1.0;
            }
        } else if constexpr (ORDER == 3) {
            switch(i) {
                case 0: return -1.0;
                case 1: return -0.447213595499958;
                case 2: return 0.447213595499958;
                case 3: return 1.0;
            }
        } else if constexpr (ORDER == 4) {
            switch(i) {
                case 0: return -1.0;
                case 1: return -0.654653670707977;
                case 2: return 0.0;
                case 3: return 0.654653670707977;
                case 4: return 1.0;
            }
        } else if constexpr (ORDER == 5) {
            switch(i) {
                case 0: return -1.0;
                case 1: return -0.765055323929465;
                case 2: return -0.285231516480645;
                case 3: return 0.285231516480645;
                case 4: return 0.765055323929465;
                case 5: return 1.0;
            }
        }
        return 0.0; // fallback for invalid i
    }
};

#endif // GLLPOINTS_H_
