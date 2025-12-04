#ifndef SRC_MODEL_MESH_IMPLEM_BUILDER_CARTESIAN_INCUDE_SEPPARAMS_H_
#define SRC_MODEL_MESH_IMPLEM_BUILDER_CARTESIAN_INCUDE_SEPPARAMS_H_

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "utils.h"

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
        this->ex_ = static_cast<ScalarType>(std::stoi(value));
      }
      else if (key == "n2")
      {
        this->ey_ = static_cast<ScalarType>(std::stoi(value));
      }
      else if (key == "n3")
      {
        this->ez_ = static_cast<ScalarType>(std::stoi(value));
      }
      else if (key == "d1")
      {
        this->hx_ = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "d2")
      {
        this->hy_ = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "d3")
      {
        this->hz_ = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "o1")
      {
        this->ox_ = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "o2")
      {
        this->oy_ = static_cast<FloatType>(std::stod(value));
      }
      else if (key == "o3")
      {
        this->oz_ = static_cast<FloatType>(std::stod(value));
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
        if (value[0] == '/' || (value.length() > 1 && value[1] == ':'))
        {
          this->data_file = value;
        }
        else
        {
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
    if (this->ex_ == utils::invalid_value<ScalarType>() ||
        this->ey_ == utils::invalid_value<ScalarType>() ||
        this->ez_ == utils::invalid_value<ScalarType>())
    {
      throw std::runtime_error(
          "Invalid SEP file: missing dimensions (n1, n2, n3)");
    }

    else if (this->ex_ <= 0 || this->ey_ <= 0 || this->ez_ <= 0)
    {
      throw std::runtime_error(
          "Invalid SEP file: invalid dimensions (n1, n2, n3)");
    }

    else if (this->hx_ == utils::invalid_value<FloatType>() ||
             this->hy_ == utils::invalid_value<FloatType>() ||
             this->hz_ == utils::invalid_value<FloatType>())
    {
      throw std::runtime_error(
          "Invalid SEP file: missing dimensions (d1, d2, d3)");
    }

    else if (this->hx_ <= 0 || this->hy_ <= 0 || this->hz_ <= 0)
    {
      throw std::runtime_error(
          "Invalid SEP file: invalid dimensions (d1, d2, d3)");
    }
  }

  /**
   * @brief Read binary data file and return as vector
   *
   * SEP format stores data in Fortran order (column-major):
   *   - n1 (x) is the fastest varying axis
   *   - n2 (y) is the middle axis
   *   - n3 (z) is the slowest varying axis
   *
   * This method reads the data and returns it with the loop order:
   *   for z in [0, n3):
   *     for y in [0, n2):
   *       for x in [0, n1):
   *         data[z * n1 * n2 + y * n1 + x]
   *
   * @tparam T Output data type (float or double)
   * @return std::vector<T> containing the model values
   * @throws std::runtime_error if file cannot be opened or read
   */
  template <typename T>
  std::vector<T> readBinaryData() const
  {
    std::ifstream file(data_file, std::ios::binary);
    if (!file.is_open())
    {
      throw std::runtime_error("Cannot open binary data file: " + data_file);
    }

    const size_t total_elements = static_cast<size_t>(ex_) * ey_ * ez_;
    std::vector<T> data(total_elements);

    // Determine if we need to convert data types
    const bool is_xdr = (data_format.find("xdr") != std::string::npos);

    if (esize == sizeof(T) && !is_xdr)
    {
      // Direct read - same size, native format
      file.read(reinterpret_cast<char*>(data.data()),
                total_elements * sizeof(T));
    }
    else if (esize == 4)
    {
      // Read as float, convert if needed
      std::vector<float> temp(total_elements);
      file.read(reinterpret_cast<char*>(temp.data()),
                total_elements * sizeof(float));

      if (is_xdr)
      {
        swapEndianness(temp.data(), total_elements);
      }

      // Convert to output type
      std::transform(temp.begin(), temp.end(), data.begin(),
                     [](float val) { return static_cast<T>(val); });
    }
    else if (esize == 8)
    {
      // Read as double, convert if needed
      std::vector<double> temp(total_elements);
      file.read(reinterpret_cast<char*>(temp.data()),
                total_elements * sizeof(double));

      if (is_xdr)
      {
        swapEndianness(temp.data(), total_elements);
      }

      // Convert to output type
      std::transform(temp.begin(), temp.end(), data.begin(),
                     [](double val) { return static_cast<T>(val); });
    }
    else
    {
      throw std::runtime_error("Unsupported esize: " + std::to_string(esize));
    }

    if (!file)
    {
      throw std::runtime_error("Error reading binary data file: " + data_file);
    }

    file.close();
    return data;
  }

  /**
   * @brief Read binary data and reorder from SEP (Fortran) to C order
   *
   * SEP stores data as [n1][n2][n3] (Fortran order, n1 fastest)
   * This returns data as [n3][n2][n1] (C order, n3 fastest... wait no)
   *
   * Actually, this reorders so the output index is:
   *   output[ix + iy * n1 + iz * n1 * n2] = sep_data[ix + iy * n1 + iz * n1 *
   * n2]
   *
   * Since SEP is already stored with x as inner loop, this matches your
   * requested order: x (inner), y (middle), z (outer)
   *
   * @tparam T Output data type
   * @return std::vector<T> with data ordered as [z][y][x] (x varies fastest)
   */
  template <typename T>
  std::vector<T> readBinaryDataXYZ() const
  {
    // SEP format already stores data with n1 (x) as fastest varying
    // So the native read order matches: for z { for y { for x { } } }
    return readBinaryData<T>();
  }

  /**
   * @brief Read binary data and reorder to ZYX order (z varies fastest)
   *
   * Reorders data so that:
   *   output[iz + iy * n3 + ix * n3 * n2]
   * corresponds to position (ix, iy, iz) in the model
   *
   * @tparam T Output data type
   * @return std::vector<T> with data ordered as [x][y][z] (z varies fastest)
   */
  template <typename T>
  std::vector<T> readBinaryDataZYX() const
  {
    std::vector<T> sep_data = readBinaryData<T>();
    std::vector<T> reordered(sep_data.size());

    const size_t nx = static_cast<size_t>(ex_);
    const size_t ny = static_cast<size_t>(ey_);
    const size_t nz = static_cast<size_t>(ez_);

    // Reorder from [z][y][x] to [x][y][z]
    for (size_t iz = 0; iz < nz; ++iz)
    {
      for (size_t iy = 0; iy < ny; ++iy)
      {
        for (size_t ix = 0; ix < nx; ++ix)
        {
          // SEP index: x + y * nx + z * nx * ny
          const size_t sep_idx = ix + iy * nx + iz * nx * ny;
          // New index: z + y * nz + x * nz * ny
          const size_t new_idx = iz + iy * nz + ix * nz * ny;
          reordered[new_idx] = sep_data[sep_idx];
        }
      }
    }

    return reordered;
  }

  /**
   * @brief Get value at specific grid position
   *
   * @param data The data vector (from readBinaryData)
   * @param ix Index along x (n1) axis
   * @param iy Index along y (n2) axis
   * @param iz Index along z (n3) axis
   * @return Value at (ix, iy, iz)
   */
  template <typename T>
  T getValue(const std::vector<T>& data, ScalarType ix, ScalarType iy,
             ScalarType iz) const
  {
    const size_t idx = static_cast<size_t>(ix) + static_cast<size_t>(iy) * ex_ +
                       static_cast<size_t>(iz) * ex_ * ey_;
    return data[idx];
  }

  /**
   * @brief Read a single XY slice at given z index
   *
   * @tparam T Output data type
   * @param iz Z index of the slice
   * @return std::vector<T> containing n1 * n2 values
   */
  template <typename T>
  std::vector<T> readSliceXY(ScalarType iz) const
  {
    if (iz < 0 || iz >= ez_)
    {
      throw std::runtime_error("Slice index out of range: " +
                               std::to_string(iz));
    }

    std::ifstream file(data_file, std::ios::binary);
    if (!file.is_open())
    {
      throw std::runtime_error("Cannot open binary data file: " + data_file);
    }

    const size_t slice_elements = static_cast<size_t>(ex_) * ey_;
    const size_t offset = static_cast<size_t>(iz) * slice_elements * esize;

    file.seekg(offset, std::ios::beg);

    std::vector<T> slice(slice_elements);
    const bool is_xdr = (data_format.find("xdr") != std::string::npos);

    if (esize == sizeof(T) && !is_xdr)
    {
      file.read(reinterpret_cast<char*>(slice.data()),
                slice_elements * sizeof(T));
    }
    else if (esize == 4)
    {
      std::vector<float> temp(slice_elements);
      file.read(reinterpret_cast<char*>(temp.data()),
                slice_elements * sizeof(float));
      if (is_xdr)
      {
        swapEndianness(temp.data(), slice_elements);
      }
      std::transform(temp.begin(), temp.end(), slice.begin(),
                     [](float val) { return static_cast<T>(val); });
    }
    else if (esize == 8)
    {
      std::vector<double> temp(slice_elements);
      file.read(reinterpret_cast<char*>(temp.data()),
                slice_elements * sizeof(double));
      if (is_xdr)
      {
        swapEndianness(temp.data(), slice_elements);
      }
      std::transform(temp.begin(), temp.end(), slice.begin(),
                     [](double val) { return static_cast<T>(val); });
    }

    if (!file)
    {
      throw std::runtime_error("Error reading slice from: " + data_file);
    }

    return slice;
  }

  void print() const
  {
    std::cout << "\n=== SEP Header Information ===\n"
              << "Dimensions: n1=" << ex_ << ", n2=" << ey_ << ", n3=" << ez_
              << "\n"
              << "Spacing:    d1=" << hx_ << ", d2=" << hy_ << ", d3=" << hz_
              << "\n"
              << "Origin:     o1=" << ox_ << ", o2=" << oy_ << ", o3=" << oz_
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
  ScalarType getN1() const { return ex_; }
  ScalarType getN2() const { return ey_; }
  ScalarType getN3() const { return ez_; }

  FloatType getD1() const { return hx_; }
  FloatType getD2() const { return hy_; }
  FloatType getD3() const { return hz_; }

  FloatType getO1() const { return ox_; }
  FloatType getO2() const { return oy_; }
  FloatType getO3() const { return oz_; }

  int getEsize() const { return esize; }
  const std::string& getDataFormat() const { return data_format; }
  const std::string& getDataFile() const { return data_file; }
  const std::string& getDataLabel() const { return data_label; }
  const std::string& getDataFiletype() const { return data_filetype; }

  int getOrder() const { return order; }
  bool isModelOnNode() const { return is_model_on_node; }
  bool isElastic() const { return is_elastic; }

  // Convenience getters for total size
  ScalarType getTotalElements() const { return ex_ * ey_ * ez_; }
  size_t getTotalBytes() const
  {
    return static_cast<size_t>(ex_) * static_cast<size_t>(ey_) *
           static_cast<size_t>(ez_) * esize;
  }

 private:
  int order;
  ScalarType ex_{utils::invalid_value<ScalarType>()},
      ey_{utils::invalid_value<ScalarType>()},
      ez_{utils::invalid_value<ScalarType>()};
  FloatType hx_{utils::invalid_value<FloatType>()},
      hy_{utils::invalid_value<FloatType>()},
      hz_{utils::invalid_value<FloatType>()};
  FloatType ox_, oy_, oz_;
  std::string data_format, data_file, data_label, data_filetype;
  int esize;
  bool is_model_on_node;
  bool is_elastic;

  /**
   * @brief Swap endianness for XDR format (big-endian to native)
   *
   * @tparam T Data type (float or double)
   * @param data Pointer to data array
   * @param count Number of elements
   */
  template <typename T>
  static void swapEndianness(T* data, size_t count)
  {
    constexpr size_t size = sizeof(T);
    auto* bytes = reinterpret_cast<char*>(data);

    for (size_t i = 0; i < count; ++i)
    {
      char* element = bytes + i * size;
      for (size_t j = 0; j < size / 2; ++j)
      {
        std::swap(element[j], element[size - 1 - j]);
      }
    }
  }

  /**
   * @brief Check if system is little-endian
   */
  static bool isLittleEndian()
  {
    const uint32_t test = 1;
    return *reinterpret_cast<const char*>(&test) == 1;
  }

  /**
   * @brief Extract directory path from a file path
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
   */
  static void trimString(std::string& str)
  {
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos)
    {
      str.clear();
      return;
    }
    str = str.substr(start);

    size_t end = str.find_last_not_of(" \t\r\n");
    if (end != std::string::npos)
    {
      str = str.substr(0, end + 1);
    }
  }
};

}  // namespace model

#endif  // SRC_MODEL_MESH_IMPLEM_BUILDER_CARTESIAN_INCUDE_SEPPARAMS_H_
