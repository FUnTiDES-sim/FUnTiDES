#ifndef FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_FIELD_SAMPLER_H_
#define FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_FIELD_SAMPLER_H_

#include <stdexcept>

#include "common_macros.h"
#include "data_type.h"

namespace solver {

/**
 * @brief Weighted sums of a nodal field at a fixed set of points, computed on the device.
 *
 * Point p takes sum_k weights(p * n + k) * field(nodes(p * n + k)), n nodes per point (typically the
 * nodes of the element containing a receiver and its interpolation weights). The nodes and weights stay
 * on the device; only the n_points results are copied back, so recording receivers at every time step
 * does not move the wavefield between device and host.
 */
class FieldSampler {
 public:
  /**
   * @param nodes Global node indices, n_points * nodes_per_point.
   * @param weights Weights, same layout.
   * @param nodes_per_point n.
   * @throws std::invalid_argument If the sizes do not match.
   */
  FieldSampler(const vectorInt& nodes, const vectorReal& weights, int nodes_per_point)
      : nodes_(nodes), weights_(weights), n_(nodes_per_point) {
    if (n_ <= 0 || nodes.extent(0) != weights.extent(0) || nodes.extent(0) % n_ != 0)
      throw std::invalid_argument("FieldSampler: nodes and weights must hold n_points * nodes_per_point values");
    n_points_ = static_cast<int>(nodes.extent(0)) / n_;
    values_ = allocateVector<vectorReal>(n_points_, "field_sampler_values");
  }

  /// Number of points.
  int numPoints() const { return n_points_; }

  /// Samples @p field; the result stays on the device (see values()).
  void sample(const vectorReal& field) const {
    auto nodes = nodes_;
    auto weights = weights_;
    auto values = values_;
    int const n = n_;
    Kokkos::parallel_for(
        "FieldSampler", Kokkos::RangePolicy<>(0, n_points_), KOKKOS_LAMBDA(const int p) {
          real_t acc = 0;
          for (int k = 0; k < n; ++k) acc += weights(p * n + k) * field(nodes(p * n + k));
          values(p) = acc;
        });
    Kokkos::fence();
  }

  /// Results of the last sample(), on the device.
  const vectorReal& values() const { return values_; }

 private:
  vectorInt nodes_;
  vectorReal weights_;
  vectorReal values_;
  int n_{0};
  int n_points_{0};
};

}  // namespace solver
#endif  // FUNTIDES_SOLVER_IMPL_COMMON_INCLUDE_FIELD_SAMPLER_H_
