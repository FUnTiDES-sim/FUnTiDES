#pragma once
#include "model_file_reader.h"

namespace model
{

/**
 * @brief Fill model Kokkos views from a .ftmd file.
 *
 * Seuls les champs présents dans le fichier sont remplis,
 * les autres sont laissés inchangés.
 */
template <typename FloatType>
void fillModelFromFile(const std::string& model_path, VECTOR_REAL_VIEW& vp_elem,
                       VECTOR_REAL_VIEW& rho_elem, VECTOR_REAL_VIEW& vs_elem,
                       VECTOR_REAL_VIEW& delta_elem,
                       VECTOR_REAL_VIEW& epsilon_elem,
                       VECTOR_REAL_VIEW& gamma_elem,
                       VECTOR_REAL_VIEW& theta_elem, VECTOR_REAL_VIEW& phi_elem)
{
  ModelFileReader reader(model_path);
  const uint64_t n = reader.nElements();

  auto fill = [&](ModelParamId id, VECTOR_REAL_VIEW& view,
                  const std::string& label) {
    if (!reader.hasParam(id)) return;
    view = allocateVector<VECTOR_REAL_VIEW>(static_cast<int>(n), label.c_str());
    const auto& buf = reader.getParam(id);
    auto host = Kokkos::create_mirror_view(view);
    for (uint64_t i = 0; i < n; ++i) host[i] = static_cast<FloatType>(buf[i]);
    Kokkos::deep_copy(view, host);
  };

  fill(ModelParamId::kVp, vp_elem, "model_vp_element");
  fill(ModelParamId::kRho, rho_elem, "model_rho_element");
  fill(ModelParamId::kVs, vs_elem, "model_vs_element");
  fill(ModelParamId::kDelta, delta_elem, "model_delta_element");
  fill(ModelParamId::kEpsilon, epsilon_elem, "model_epsilon_element");
  fill(ModelParamId::kGamma, gamma_elem, "model_gamma_element");
  fill(ModelParamId::kTheta, theta_elem, "model_theta_element");
  fill(ModelParamId::kPhi, phi_elem, "model_phi_element");
}

}  // namespace model