#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace model
{
template <typename FloatType, typename ScalarType>
class SepParams
{
 public:
  SepParams(const std::string& filename)
      : is_model_on_node(true), is_elastic(false)
  {
    readSepFile(filename);
  }

  void readSepFile(const std::string& filename)
  {
    std::ifstream file(filename);

    if (!file.is_open())
    {
      throw std::runtime_error("Cannot open SEP file: " + filename);
    }

    std::string line;
    std::string directory = getDirectory(filename);

    while (std::getline(file, line))
    {
      // Skip empty lines and comments
      if (line.empty() || line[0] == '#')
      {
        continue;
      }

      // Trim leading whitespace
      size_t start = line.find_first_not_of(" \t");
      if (start == std::string::npos)
      {
        continue;
      }

      line = line.substr(start);

      // Split by '='
      size_t eq_pos = line.find('=');
      if (eq_pos == std::string::npos)
      {
        continue;
      }

      std::string key = line.substr(0, eq_pos);
      std::string value = line.substr(eq_pos + 1);

      // Trim whitespace from key and value
      trimString(key);
      trimString(value);

      // Parse parameters
      if (key == "n1")
      {
        this->ex = static_cast<ScalarType>(std::stoi(value));
      }
      else if (key == "n2")
      {
        this->ey = static_cast<ScalarType>(std::stoi(value));
      }
      else if (key == "n3")
      {
        this->ez = static_cast<ScalarType>(std::stoi(value));
      }
      else if (key == "d1")
      {
        this->hx = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "d2")
      {
        this->hy = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "d3")
      {
        this->hz = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "o1")
      {
        this->ox = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "o2")
      {
        this->oy = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "o3")
      {
        this->oz = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "esize")
      {
        this->esize = std::stoi(value);
      }
      else if (key == "data_format")
      {
        this->data_format = value;
      }
      else if (key == "in")
      {
        // Handle data file path - may be relative
        if (value[0] == '/' || (value.length() > 1 && value[1] == ':'))
        {
          // Absolute path
          this->data_file = value;
        }
        else
        {
          // Relative path - resolve relative to header directory
          this->data_file = directory + "/" + value;
        }
      }
      else if (key == "data_label")
      {
        this->data_label = value;
      }
      else if (key == "data_filetype")
      {
        this->data_filetype = value;
      }
      else if (key == "order")
      {
        this->order = std::stoi(value);
      }
      else if (key == "model_on_node")
      {
        this->is_model_on_node = (value == "true" || value == "1");
      }
      else if (key == "elastic")
      {
        this->is_elastic = (value == "true" || value == "1");
      }
    }

    file.close();

    // Validate critical parameters
    if (this->ex == 0 || this->ey == 0 || this->ez == 0)
    {
      throw std::runtime_error(
          "Invalid SEP file: missing or invalid dimensions (n1, n2, n3)");
    }
  }

  void print() const
  {
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
  size_t getTotalBytes() const
  {
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
  static std::string getDirectory(const std::string& filepath)
  {
    size_t pos = filepath.find_last_of("/\\");
    if (pos == std::string::npos)
    {
      return ".";
    }
    if (pos == 0)
    {
      return "/";
    }
    return filepath.substr(0, pos);
  }

  /**
   * @brief Trim leading and trailing whitespace from a string in-place
   *
   * @param str String to trim
   */
  static void trimString(std::string& str)
  {
    // Trim leading whitespace
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos)
    {
      str.clear();
      return;
    }
    str = str.substr(start);

    // Trim trailing whitespace
    size_t end = str.find_last_not_of(" \t\r\n");
    if (end != std::string::npos)
    {
      str = str.substr(0, end + 1);
    }
  }
};
}  // namespace model
