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
 * @brief Parses a text model file and exposes one value array per property.
 *
 * File format: one value per line, one section per property.
 * @code
 * Model Vp element
 * <number of values>
 * value 1
 * ...
 * value N
 *
 * Model Rho element
 * <number of values>
 * value 1
 * ...
 * value N
 * @endcode
 *
 * The support keyword of a section is `element` or `node`, and all sections of
 * a file must use the same one. Every section must hold the same number of
 * values, and a property may appear only once. Lines that do not start with
 * "Model" outside a section are ignored. Property names are stored as written
 * in the file (for example Vp, Vs, Rho, Qp, Qs, Delta, Epsilon, Gamma, Theta,
 * Phi).
 * @todo VERIFY: are property names case-sensitive for the consumers, and are
 * the units of each property (for example Theta and Phi) defined by the file
 * or by the consumer?
 */
class CartesianModelFileReader {
 public:
  /**
   * @brief Constructs the reader and parses the whole file.
   * @param path Path to the model text file.
   * @throws std::runtime_error if the file cannot be opened or is malformed.
   */
  explicit CartesianModelFileReader(const std::string& path) { parse(path); }

  /**
   * @brief Tells whether a property is present in the file.
   * @param prop Property name, as written in the section header.
   * @return true if a section for prop was parsed.
   */
  bool has(const std::string& prop) const { return data_.find(prop) != data_.end(); }

  /**
   * @brief Returns the values of a property, in file order.
   * @param prop Property name, as written in the section header.
   * @return Array of count() values.
   * @throws std::runtime_error if the property is absent.
   */
  const std::vector<double>& get(const std::string& prop) const {
    auto it = data_.find(prop);
    if (it == data_.end())
      throw std::runtime_error("CartesianModelFileReader: property '" + prop + "' not found in file");
    return it->second;
  }

  /**
   * @brief Tells whether the values are attached to nodes or to elements.
   * @return true if the support is `node`, false if it is `element`.
   */
  bool onNodes() const { return on_nodes_; }

  /**
   * @brief Returns the number of values of each property.
   */
  size_t count() const { return count_; }

 private:
  std::map<std::string, std::vector<double>> data_;  ///< Values by property name.
  bool on_nodes_{false};                             ///< Support shared by all sections.
  size_t count_{0};                                  ///< Values per property.
  bool support_set_{false};                          ///< True once the first section fixed the support.

  void parse(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) throw std::runtime_error("CartesianModelFileReader: cannot open file '" + path + "'");

    for (std::string line; std::getline(file, line);) {
      if (line.empty()) continue;

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
        // Blank lines between the header and the count are allowed.
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
