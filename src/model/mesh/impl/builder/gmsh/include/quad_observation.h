#ifndef FUNTIDES_QUAD_OBSERVATION_H_
#define FUNTIDES_QUAD_OBSERVATION_H_
#include <cmath>
#include <stdexcept>
#include <vector>

#include "data_type.h"
#include "gllpoints.h"
#include "model_unstruct.h"

namespace model {
struct QuadObservation {
  int element;
  double x, z;
  std::vector<int> nodes;
  std::vector<double> weights;
};

/// Locate an observation point by inversion of the bilinear Q4 geometry, then
/// evaluate the SEM Lagrange basis. Never snaps a receiver to a nearby node.
inline QuadObservation locateQuadObservation(const ModelUnstruct<float, int>& mesh, double x, double z) {
  int q = mesh.getOrder(), n = q + 1;
  for (int e = 0; e < mesh.getNumberOfElements(); ++e) {
    int ids[4] = {mesh.globalNodeIndex(e, 0, 0, 0), mesh.globalNodeIndex(e, q, 0, 0), mesh.globalNodeIndex(e, q, 0, q),
                  mesh.globalNodeIndex(e, 0, 0, q)};
    double X[4], Z[4], minx = 1e300, maxx = -1e300, minz = 1e300, maxz = -1e300;
    for (int k = 0; k < 4; ++k) {
      X[k] = mesh.nodeCoord(ids[k], 0);
      Z[k] = mesh.nodeCoord(ids[k], 2);
      minx = std::min(minx, X[k]);
      maxx = std::max(maxx, X[k]);
      minz = std::min(minz, Z[k]);
      maxz = std::max(maxz, Z[k]);
    }
    if (x < minx - 1e-7 || x > maxx + 1e-7 || z < minz - 1e-7 || z > maxz + 1e-7) continue;
    double xi = 0, eta = 0;
    bool converged = false;
    for (int iteration = 0; iteration < 30; ++iteration) {
      double w[4] = {(1 - xi) * (1 - eta) / 4, (1 + xi) * (1 - eta) / 4, (1 + xi) * (1 + eta) / 4,
                     (1 - xi) * (1 + eta) / 4};
      double wx[4] = {-(1 - eta) / 4, (1 - eta) / 4, (1 + eta) / 4, -(1 + eta) / 4};
      double wz[4] = {-(1 - xi) / 4, -(1 + xi) / 4, (1 + xi) / 4, (1 - xi) / 4};
      double a = 0, b = 0, c = 0, d = 0, rx = -x, rz = -z;
      for (int k = 0; k < 4; ++k) {
        rx += w[k] * X[k];
        rz += w[k] * Z[k];
        a += wx[k] * X[k];
        b += wz[k] * X[k];
        c += wx[k] * Z[k];
        d += wz[k] * Z[k];
      }
      if (std::hypot(rx, rz) < 1e-7) {
        converged = true;
        break;
      }
      double det = a * d - b * c;
      if (!(det > 0)) break;
      xi -= (d * rx - b * rz) / det;
      eta -= (-c * rx + a * rz) / det;
    }
    if (!converged || std::abs(xi) > 1 + 1e-7 || std::abs(eta) > 1 + 1e-7) continue;
    QuadObservation result{e, x, z, {}, {}};
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        double wi = 1, wj = 1;
        for (int k = 0; k < n; ++k) {
          if (k != i) wi *= (xi - GLLPoints::get(q, k)) / (GLLPoints::get(q, i) - GLLPoints::get(q, k));
          if (k != j) wj *= (eta - GLLPoints::get(q, k)) / (GLLPoints::get(q, j) - GLLPoints::get(q, k));
        }
        result.nodes.push_back(mesh.globalNodeIndex(e, i, 0, j));
        result.weights.push_back(wi * wj);
      }
    return result;
  }
  throw std::runtime_error("Observation outside the quadrilateral mesh");
}
inline double observeQuad(const QuadObservation& observation, const std::vector<float>& field) {
  double value = 0;
  for (size_t i = 0; i < observation.nodes.size(); ++i)
    value += observation.weights[i] * field.at(observation.nodes[i]);
  return value;
}
}  // namespace model
#endif
