#pragma once
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "data_type.h"

namespace model
{

/**
 * @brief Parameter identifiers — must match Python PARAM_IDS dict
 */
enum class ModelParamId : uint32_t
{
  kVp = 0,
  kRho = 1,
  kVs = 2,
  kDelta = 3,
  kEpsilon = 4,
  kGamma = 5,
  kTheta = 6,
  kPhi = 7,
};

/**
 * @brief Reader for .ftmd binary model files.
 *
 * Layout:
 *   magic[4] | version(u32) | n_elem(u64) | n_params(u32)
 *   param_ids[n_params x u32]
 *   data[n_params x n_elem x f64]  (param-major)
 */
class ModelFileReader
{
 public:
  static constexpr const char* kMagic = "FTMD";

  explicit ModelFileReader(const std::string& path) { open(path); }

  uint64_t nElements() const { return n_elem_; }

  bool hasParam(ModelParamId id) const { return data_.count(id) > 0; }

  const std::vector<double>& getParam(ModelParamId id) const
  {
    auto it = data_.find(id);
    if (it == data_.end())
      throw std::runtime_error("[ModelFileReader] Parameter " +
                               std::to_string(static_cast<uint32_t>(id)) +
                               " not found in model file.");
    return it->second;
  }

 private:
  uint64_t n_elem_{0};
  std::unordered_map<ModelParamId, std::vector<double>> data_;

  void open(const std::string& path)
  {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open())
      throw std::runtime_error("[ModelFileReader] Cannot open: " + path);

    // Magic
    char magic[4];
    f.read(magic, 4);
    if (std::string(magic, 4) != kMagic)
      throw std::runtime_error("[ModelFileReader] Bad magic in: " + path);

    // Header
    uint32_t version, n_params;
    f.read(reinterpret_cast<char*>(&version), sizeof(version));
    f.read(reinterpret_cast<char*>(&n_elem_), sizeof(n_elem_));
    f.read(reinterpret_cast<char*>(&n_params), sizeof(n_params));

    // Param IDs
    std::vector<ModelParamId> ids(n_params);
    for (auto& id : ids)
    {
      uint32_t raw;
      f.read(reinterpret_cast<char*>(&raw), sizeof(raw));
      id = static_cast<ModelParamId>(raw);
    }

    // Data (param-major)
    for (auto id : ids)
    {
      std::vector<double> buf(n_elem_);
      f.read(reinterpret_cast<char*>(buf.data()),
             static_cast<std::streamsize>(n_elem_ * sizeof(double)));
      if (!f)
        throw std::runtime_error("[ModelFileReader] Truncated data for param " +
                                 std::to_string(static_cast<uint32_t>(id)));
      data_[id] = std::move(buf);
    }
  }
};

}  // namespace model