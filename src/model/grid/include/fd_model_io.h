/**
 * @file fdtd_model_io.h
 * @brief Model file I/O operations
 */
#ifndef SRC_MODEL_GRID_INCLUDE_FDTD_MODEL_IO_H_
#define SRC_MODEL_GRID_INCLUDE_FDTD_MODEL_IO_H_

#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "fd_grid_geometry.h"

namespace model
{
namespace fdgrid
{

/**
 * @brief Type traits for detecting container capabilities
 */
namespace detail
{

// Check if type has resize() method (std::vector)
template <typename T, typename = void>
struct has_resize : std::false_type
{
};

template <typename T>
struct has_resize<T, std::void_t<decltype(std::declval<T>().resize(size_t{}))>>
    : std::true_type
{
};

// Check if type has data() method
template <typename T, typename = void>
struct has_data : std::false_type
{
};

template <typename T>
struct has_data<T, std::void_t<decltype(std::declval<T>().data())>>
    : std::true_type
{
};

}  // namespace detail

/**
 * @brief Model file I/O utilities
 *
 * Handles reading and writing velocity models from/to binary files.
 * Supports both std::vector and Kokkos::View containers.
 *
 * **File Format:**
 * - Header: 6 values (3 ints: nx,ny,nz; 3 floats: dx,dy,dz)
 * - Data: nx*ny*nz floats (velocity values)
 */
class ModelIO
{
 public:
  /**
   * @brief Read grid geometry from model file
   * @param file_path Path to binary model file
   * @return GridGeometry object with dimensions and spacing
   * @throws std::runtime_error if file cannot be read
   */
  static GridGeometry ReadGeometry(const std::string& file_path)
  {
    std::ifstream infile(file_path, std::ios::in | std::ios::binary);
    if (!infile)
    {
      throw std::runtime_error("Cannot open model file: " + file_path);
    }

    int nx, ny, nz;
    float dx, dy, dz;

    infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nz), sizeof(int));
    infile.read(reinterpret_cast<char*>(&dx), sizeof(float));
    infile.read(reinterpret_cast<char*>(&dy), sizeof(float));
    infile.read(reinterpret_cast<char*>(&dz), sizeof(float));

    if (!infile)
    {
      throw std::runtime_error("Corrupted model file header: " + file_path);
    }

    printf("Loaded geometry from %s: %dx%dx%d, spacing %.3fx%.3fx%.3f m\n",
           file_path.c_str(), nx, ny, nz, dx, dy, dz);

    return GridGeometry(nx, ny, nz, dx, dy, dz);
  }

  /**
   * @brief Read velocity model data from file
   *
   * Works with both std::vector (has resize()) and Kokkos::View (no resize()).
   * For Kokkos Views, the view must be pre-allocated with correct size.
   *
   * @param file_path Path to binary model file
   * @param geom Grid geometry (must match file)
   * @param data Output container (std::vector or Kokkos::View)
   * @throws std::runtime_error if read fails or size mismatch
   */
  template <typename VectorType>
  static void ReadVelocityModel(const std::string& file_path,
                                const GridGeometry& geom, VectorType& data)
  {
    std::ifstream infile(file_path, std::ios::in | std::ios::binary);
    if (!infile)
    {
      throw std::runtime_error("Cannot open model file: " + file_path);
    }

    // Skip header (6 values: 3 ints + 3 floats)
    infile.seekg(3 * sizeof(int) + 3 * sizeof(float));

    size_t expected_size = geom.TotalPoints();

    // Resize if container supports it (std::vector)
    if constexpr (detail::has_resize<VectorType>::value)
    {
      data.resize(expected_size);
    }
    else
    {
      // For Kokkos::View, check size matches
      if (data.extent(0) != expected_size)
      {
        throw std::runtime_error("Kokkos::View size mismatch: expected " +
                                 std::to_string(expected_size) + ", got " +
                                 std::to_string(data.extent(0)) +
                                 ". Pre-allocate view with correct size.");
      }
    }

    // Read into temporary buffer for Kokkos (host -> device transfer)
    std::vector<float> temp_buffer(expected_size);
    infile.read(reinterpret_cast<char*>(temp_buffer.data()),
                expected_size * sizeof(float));

    if (!infile)
    {
      throw std::runtime_error("Failed to read velocity data: " + file_path);
    }

    // Copy to output (handles both std::vector and Kokkos::View)
    CopyData(temp_buffer, data, expected_size);

    printf("Loaded %zu velocity values from %s\n", expected_size,
           file_path.c_str());
  }

  /**
   * @brief Write velocity model to file
   * @param file_path Output file path
   * @param geom Grid geometry
   * @param data Velocity model data (std::vector or Kokkos::View)
   * @throws std::runtime_error if write fails
   */
  template <typename VectorType>
  static void WriteVelocityModel(const std::string& file_path,
                                 const GridGeometry& geom,
                                 const VectorType& data)
  {
    std::ofstream outfile(file_path, std::ios::out | std::ios::binary);
    if (!outfile)
    {
      throw std::runtime_error("Cannot create model file: " + file_path);
    }

    // Write header
    int nx = geom.nx(), ny = geom.ny(), nz = geom.nz();
    float dx = geom.dx(), dy = geom.dy(), dz = geom.dz();
    outfile.write(reinterpret_cast<const char*>(&nx), sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&ny), sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&nz), sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&dx), sizeof(float));
    outfile.write(reinterpret_cast<const char*>(&dy), sizeof(float));
    outfile.write(reinterpret_cast<const char*>(&dz), sizeof(float));

    // Copy to temporary buffer if needed (Kokkos device -> host)
    std::vector<float> temp_buffer;
    const float* write_ptr = GetDataPointer(data, temp_buffer);

    // Write data
    outfile.write(reinterpret_cast<const char*>(write_ptr),
                  data.extent(0) * sizeof(float));

    if (!outfile)
    {
      throw std::runtime_error("Failed to write model file: " + file_path);
    }

    printf("Wrote velocity model to %s\n", file_path.c_str());
  }

 private:
  /**
   * @brief Copy data from source to destination (handles different container
   * types)
   */
  template <typename SourceType, typename DestType>
  static void CopyData(const SourceType& src, DestType& dest, size_t size)
  {
    // Direct copy for std::vector or containers with direct data access
    if constexpr (detail::has_data<DestType>::value)
    {
      std::copy(src.begin(), src.begin() + size, dest.data());
    }
    else
    {
      // Element-wise copy for Kokkos::View (triggers host->device transfer)
      for (size_t i = 0; i < size; ++i)
      {
        dest(i) = src[i];
      }
    }
  }

  /**
   * @brief Get raw data pointer, using temporary buffer for Kokkos::View if
   * needed
   */
  template <typename VectorType>
  static const float* GetDataPointer(const VectorType& data,
                                     std::vector<float>& temp_buffer)
  {
    if constexpr (detail::has_data<VectorType>::value)
    {
      // std::vector or similar - direct access
      return data.data();
    }
    else
    {
      // Kokkos::View - copy to host buffer first
      temp_buffer.resize(data.extent(0));
      for (size_t i = 0; i < data.extent(0); ++i)
      {
        temp_buffer[i] = data(i);
      }
      return temp_buffer.data();
    }
  }
};

}  // namespace fdgrid
}  // namespace model

#endif  // SRC_MODEL_GRID_INCLUDE_FDTD_MODEL_IO_H_
