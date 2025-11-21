#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace model {
template <typename FloatType, typename ScalarType>
struct SepParams {
 public:
  SepParams()
      : order(1),
        ex(0),
        ey(0),
        ez(0),
        hx(1.0),
        hy(1.0),
        hz(1.0),
        ox(0.0),
        oy(0.0),
        oz(0.0),
        esize(4),
        is_model_on_node(true),
        is_elastic(false) {}

  static SepParams readSepFile(const std::string& filename) {
    SepParams params;
    std::ifstream file(filename);

    if (!file.is_open()) {
      throw std::runtime_error("Cannot open SEP file: " + filename);
    }

    std::string line;
    std::string directory = getDirectory(filename);

    while (std::getline(file, line)) {
      // Skip empty lines and comments
      if (line.empty() || line[0] == '#') {
        continue;
      }

      // Trim leading whitespace
      size_t start = line.find_first_not_of(" \t");
      if (start == std::string::npos) {
        continue;
      }

      line = line.substr(start);

      // Split by '='
      size_t eq_pos = line.find('=');
      if (eq_pos == std::string::npos) {
        continue;
      }

      std::string key = line.substr(0, eq_pos);
      std::string value = line.substr(eq_pos + 1);

      // Trim whitespace from key and value
      trimString(key);
      trimString(value);

      // Parse parameters
      if (key == "n1") {
        params.ex = static_cast<ScalarType>(std::stoi(value));
      } else if (key == "n2") {
        params.ey = static_cast<ScalarType>(std::stoi(value));
      } else if (key == "n3") {
        params.ez = static_cast<ScalarType>(std::stoi(value));
      } else if (key == "d1") {
        params.hx = static_cast<FloatType>(std::stod(value));
      } else if (key == "d2") {
        params.hy = static_cast<FloatType>(std::stod(value));
      } else if (key == "d3") {
        params.hz = static_cast<FloatType>(std::stod(value));
      } else if (key == "o1") {
        params.ox = static_cast<FloatType>(std::stod(value));
      } else if (key == "o2") {
        params.oy = static_cast<FloatType>(std::stod(value));
      } else if (key == "o3") {
        params.oz = static_cast<FloatType>(std::stod(value));
      } else if (key == "esize") {
        params.esize = std::stoi(value);
      } else if (key == "data_format") {
        params.data_format = value;
      } else if (key == "in") {
        // Handle data file path - may be relative
        if (value[0] == '/' || (value.length() > 1 && value[1] == ':')) {
          // Absolute path
          params.data_file = value;
        } else {
          // Relative path - resolve relative to header directory
          params.data_file = directory + "/" + value;
        }
      } else if (key == "data_label") {
        params.data_label = value;
      } else if (key == "data_filetype") {
        params.data_filetype = value;
      } else if (key == "order") {
        params.order = std::stoi(value);
      } else if (key == "model_on_node") {
        params.is_model_on_node = (value == "true" || value == "1");
      } else if (key == "elastic") {
        params.is_elastic = (value == "true" || value == "1");
      }
    }

    file.close();

    // Validate critical parameters
    if (params.ex == 0 || params.ey == 0 || params.ez == 0) {
      throw std::runtime_error(
          "Invalid SEP file: missing or invalid dimensions (n1, n2, n3)");
    }

    return params;
  }

  void print() const {
    std::cout << "\n=== SEP Header Information ===\n"
              << "Dimensions: n1=" << ex << ", n2=" << ey << ", n3=" << ez
              << "\n"
              << "Spacing:    d1=" << hx << ", d2=" << hy << ", d3=" << hz
              << "\n"
              << "Origin:     o1=" << ox << ", o2=" << oy << ", o3=" << oz
              << "\n"
              << "Data format: " << data_format << " (esize=" << esize
              << " bytes)\n"
              << "Data file:  " << data_file << "\n"
              << "Data label: " << data_label << "\n"
              << "File type:  " << data_filetype << "\n"
              << "FE Order:   " << order << "\n"
              << "Model on node: " << (is_model_on_node ? "true" : "false")
              << "\n"
              << "Elastic: " << (is_elastic ? "true" : "false") << "\n"
              << "==============================\n\n";
  }

  // Getter methods
  ScalarType getN1() const { return ex; }
  ScalarType getN2() const { return ey; }
  ScalarType getN3() const { return ez; }

  FloatType getD1() const { return hx; }
  FloatType getD2() const { return hy; }
  FloatType getD3() const { return hz; }

  FloatType getO1() const { return ox; }
  FloatType getO2() const { return oy; }
  FloatType getO3() const { return oz; }

  int getEsize() const { return esize; }
  const std::string& getDataFormat() const { return data_format; }
  const std::string& getDataFile() const { return data_file; }
  const std::string& getDataLabel() const { return data_label; }
  const std::string& getDataFiletype() const { return data_filetype; }

  int getOrder() const { return order; }
  bool isModelOnNode() const { return is_model_on_node; }
  bool isElastic() const { return is_elastic; }

  // Convenience getters for total size
  ScalarType getTotalElements() const { return ex * ey * ez; }
  size_t getTotalBytes() const {
    return static_cast<size_t>(ex) * static_cast<size_t>(ey) *
           static_cast<size_t>(ez) * esize;
  }

 private:
  int order;
  ScalarType ex, ey, ez;
  FloatType hx, hy, hz;
  FloatType ox, oy, oz;
  std::string data_format, data_file, data_label, data_filetype;
  int esize;
  bool is_model_on_node;
  bool is_elastic;

  /**
   * @brief Extract directory path from a file path
   *
   * @param filepath Full or relative file path
   *
   * @return Directory component of the path
   *   - Returns "." for files in current directory
   *   - Returns "/" for root directory files
   *   - Returns the directory path otherwise
   *
   * @details
   * Used internally to resolve relative paths of data files relative to
   * the directory containing the header file. Handles both Unix-style (/)
   * and Windows-style (\\) path separators.
   */
  static std::string getDirectory(const std::string& filepath) {
    size_t pos = filepath.find_last_of("/\\");
    if (pos == std::string::npos) {
      return ".";
    }
    if (pos == 0) {
      return "/";
    }
    return filepath.substr(0, pos);
  }

  /**
   * @brief Trim leading and trailing whitespace from a string in-place
   *
   * @param str String to trim
   */
  static void trimString(std::string& str) {
    // Trim leading whitespace
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
      str.clear();
      return;
    }
    str = str.substr(start);

    // Trim trailing whitespace
    size_t end = str.find_last_not_of(" \t\r\n");
    if (end != std::string::npos) {
      str = str.substr(0, end + 1);
    }
  }
};
}  // namespace model
