#ifndef FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_MODEL_FILE_READER_H_
#define FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_CARTESIAN_MODEL_FILE_READER_H_

#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace model {

/**
 * @brief Parses a text model file and exposes per-property value arrays.
 *
 * File format (one value per line):
 * @code
 * Model Vp element
 * <Nber of elements>
 * Value 1
 * Value 2
 * ...
 * Value N
 *
 * Model Rho element
 * <Nber of elements>
 * Value 1
 * Value 2
 * ...
 * Value N
 * @endcode
 *
 * All sections in a file must share the same support (`element` or `node`).
 * Supported property names: Vp, Vs, Rho, Qp, Qs, Delta, Epsilon, Gamma,
 * Theta, Phi.
 */
class CartesianModelFileReader {
 public:
  /**
   * @brief Constructs the reader and immediately parses the file.
   * @param path Path to the model text file.
   * @throws std::runtime_error if the file cannot be opened or is malformed.
   */
  explicit CartesianModelFileReader(const std::string& path) { parse(path); }

  /**
   * @brief Returns whether a property was found in the file.
   * @param prop Property name (e.g. "Vp", "Vs", "Rho").
   */
  bool has(const std::string& prop) const { return data_.find(prop) != data_.end(); }

  /**
   * @brief Returns the parsed values for a given property.
   * @param prop Property name.
   * @throws std::runtime_error if the property is absent.
   */
  const std::vector<double>& get(const std::string& prop) const {
    auto it = data_.find(prop);
    if (it == data_.end())
      throw std::runtime_error("CartesianModelFileReader: property '" + prop + "' not found in file");
    return it->second;
  }

  /**
   * @brief Returns true if the support is `node`, false if `element`.
   */
  bool onNodes() const { return on_nodes_; }

  /**
   * @brief Returns the number of values per property.
   */
  size_t count() const { return count_; }

 private:
  std::map<std::string, std::vector<double>> data_;
  bool on_nodes_{false};
  size_t count_{0};
  bool support_set_{false};

  void parse(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) throw std::runtime_error("CartesianModelFileReader: cannot open file '" + path + "'");

    for (std::string line; std::getline(file, line);) {
      if (line.empty()) continue;

      // Detect section header: "Model <Prop> <support>"
      std::stringstream header(line);
      std::string keyword;
      header >> keyword;
      if (keyword != "Model") continue;

      std::string prop, support;
      if (!(header >> prop >> support))
        throw std::runtime_error("CartesianModelFileReader: malformed section header: '" + line + "'");

      bool section_on_nodes = (support == "node");
      if (support != "node" && support != "element")
        throw std::runtime_error("CartesianModelFileReader: unknown support '" + support +
                                 "' (expected 'element' or 'node')");

      if (!support_set_) {
        on_nodes_ = section_on_nodes;
        support_set_ = true;
      } else if (section_on_nodes != on_nodes_) {
        throw std::runtime_error("CartesianModelFileReader: mixed supports in file");
      }

      std::string count_line;
      for (; std::getline(file, count_line) && count_line.empty();) {
        // Skip empty lines between header and count
      }

      size_t section_count = 0;
      try {
        section_count = static_cast<size_t>(std::stoul(count_line));
      } catch (...) {
        throw std::runtime_error("CartesianModelFileReader: expected integer count after 'Model " + prop + "', got: '" +
                                 count_line + "'");
      }

      if (count_ == 0)
        count_ = section_count;
      else if (section_count != count_)
        throw std::runtime_error("CartesianModelFileReader: property count mismatch for '" + prop + "'");

      std::vector<double> values;
      values.reserve(section_count);

      for (std::string val_line; values.size() < section_count;) {
        if (!std::getline(file, val_line))
          throw std::runtime_error("CartesianModelFileReader: unexpected end of file while reading '" + prop + "'");

        if (val_line.empty()) continue;

        try {
          values.push_back(std::stod(val_line));
        } catch (...) {
          throw std::runtime_error("CartesianModelFileReader: non-numeric value '" + val_line + "' in '" + prop + "'");
        }
      }

      if (data_.count(prop)) throw std::runtime_error("CartesianModelFileReader: duplicate property '" + prop + "'");
      data_[prop] = std::move(values);
    }

    if (data_.empty()) throw std::runtime_error("CartesianModelFileReader: no valid section found");
  };
};

}  // namespace model
#endif  // FUNTIDES_MODEL_MESH_IMPL_BUILDER_CARTESIAN_INCLUDE_MODEL_FILE_READER_H_
