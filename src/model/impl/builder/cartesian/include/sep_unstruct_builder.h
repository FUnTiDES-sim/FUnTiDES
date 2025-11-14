/**
 * @file sep_unstruct_builder.h
 * @brief SEP file-based builder for unstructured spectral element meshes
 *
 * This file contains the SepUnstructBuilder class which constructs unstructured
 * spectral element meshes from SEP (Stanford Exploration Project) format files.
 * The builder supports both acoustic (velocity-only) and elastic
 * (velocity+density) models with spectral orders from 1 to 5.
 */

#ifndef SRC_MODEL_SEP_INCLUDE_SEP_UNSTRUCT_BUILDER_H_
#define SRC_MODEL_SEP_INCLUDE_SEP_UNSTRUCT_BUILDER_H_

#include <builder.h>
#include <data_type.h>
#include <model_unstruct.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace model {

/**
 * @class SepUnstructBuilder
 * @brief Builder class for creating unstructured meshes from SEP files
 *
 * This class reads velocity and optional density models from SEP format files
 * and constructs a complete unstructured spectral element mesh. The SEP format
 * consists of an ASCII header file containing metadata (dimensions, sampling
 * intervals) and a binary data file containing the actual physical property
 * values.
 *
 * The builder performs the following operations:
 * 1. Parses SEP header files to extract grid geometry
 * 2. Loads binary data for velocity and optionally density
 * 3. Constructs spectral element mesh with GLL nodes
 * 4. Maps physical properties to mesh nodes or elements
 *
 * @tparam FloatType Floating point type for physical properties (e.g., double,
 * float)
 * @tparam ScalarType Scalar integer type for dimensions (e.g., int, size_t)
 *
 * @note SEP files store data in element-centered format (one value per grid
 * cell). This builder handles conversion to node-centered format when required.
 *
 * @example
 * // Create acoustic model (velocity only)
 * auto builder = model::SepUnstructBuilder<double, int>(
 *     "vp_model.H", 4, true);
 * auto mesh = builder.getModel();
 *
 * @example
 * // Create elastic model (velocity + density)
 * auto builder = model::SepUnstructBuilder<double, int>(
 *     "vp_model.H", "rho_model.H", 3, false);
 * auto mesh = builder.getModel();
 */
template <typename FloatType, typename ScalarType>
class SepUnstructBuilder : public ModelBuilderBase<FloatType, ScalarType> {
 public:
  using ModelBuilderBase<FloatType, ScalarType>::MAX_ORDER;

  /**
   * @brief Constructs builder with velocity and optional density models
   *
   * @param vp_sep_file Path to SEP header file containing velocity model
   * @param rho_sep_file Path to SEP header file containing density model
   *                     (empty string for acoustic models)
   * @param order Spectral element order (1-5)
   * @param is_model_on_nodes If true, stores model on nodes; if false, on
   * elements
   *
   * @throws std::runtime_error if order is outside valid range [1, MAX_ORDER]
   * @throws std::runtime_error if SEP files cannot be opened
   * @throws std::runtime_error if density dimensions don't match velocity
   * @throws std::runtime_error if required SEP parameters are missing
   *
   * @note Constructor performs all initialization including mesh generation
   *       and data loading. The model is ready to use after construction.
   */
  SepUnstructBuilder(const std::string& vp_sep_file,
                     const std::string& rho_sep_file, int order,
                     bool is_model_on_nodes)
      : vp_sep_file_(vp_sep_file),
        rho_sep_file_(rho_sep_file),
        order_(order),
        is_model_on_nodes_(is_model_on_nodes),
        has_density_(!rho_sep_file.empty()) {
    if (order_ < 1 || order_ > MAX_ORDER) {
      throw std::runtime_error("Order must be between 1 and " +
                               std::to_string(MAX_ORDER));
    }

    ParseSepHeader(vp_sep_file_);
    ExtractModelParameters();
    InitGlobalNodeList();
    InitNodesCoords();
    LoadVelocity();
    if (has_density_) {
      ValidateDensityGeometry();
      LoadDensity();
    } else {
      InitDefaultDensity();
    }
  }

  /**
   * @brief Constructs acoustic model builder (velocity only)
   *
   * This is a convenience constructor for acoustic models that don't require
   * density information. Default density of 1.0 is assigned to all nodes or
   * elements.
   *
   * @param vp_sep_file Path to SEP header file containing velocity model
   * @param order Spectral element order (1-5)
   * @param is_model_on_nodes If true, stores model on nodes; if false, on
   * elements
   *
   * @throws std::runtime_error if order is outside valid range [1, MAX_ORDER]
   * @throws std::runtime_error if SEP file cannot be opened
   * @throws std::runtime_error if required SEP parameters are missing
   */
  SepUnstructBuilder(const std::string& vp_sep_file, int order,
                     bool is_model_on_nodes)
      : SepUnstructBuilder(vp_sep_file, "", order, is_model_on_nodes) {}

  /**
   * @brief Destructor
   */
  ~SepUnstructBuilder() = default;

  /**
   * @brief Constructs and returns the complete model
   *
   * Creates a ModelUnstruct instance with all mesh data, node coordinates,
   * and physical properties (velocity and density). The returned model is
   * ready for use in simulations.
   *
   * @return Shared pointer to constructed ModelApi instance
   *
   * @note This method can be called multiple times and will return a new
   *       model instance each time with the same data.
   */
  std::shared_ptr<model::ModelApi<FloatType, ScalarType>> getModel()
      const override {
    model::ModelUnstructData<FloatType, ScalarType> model_data;

    model_data.order_ = order_;
    model_data.n_element_ = ex_ * ey_ * ez_;
    model_data.n_node_ =
        (ex_ * order_ + 1) * (ey_ * order_ + 1) * (ez_ * order_ + 1);
    model_data.lx_ = lx_;
    model_data.ly_ = ly_;
    model_data.lz_ = lz_;

    model_data.global_node_index_ = global_node_index_;
    model_data.nodes_coords_x_ = nodes_coords_x_;
    model_data.nodes_coords_y_ = nodes_coords_y_;
    model_data.nodes_coords_z_ = nodes_coords_z_;

    model_data.isModelOnNodes_ = is_model_on_nodes_;

    model_data.model_vp_node_ = model_vp_node_;
    model_data.model_rho_node_ = model_rho_node_;
    model_data.model_vp_element_ = model_vp_element_;
    model_data.model_rho_element_ = model_rho_element_;

    return std::make_shared<model::ModelUnstruct<FloatType, ScalarType>>(
        model_data);
  }

  /**
   * @brief Gets number of elements in x direction
   * @return Number of elements along x-axis
   */
  ScalarType GetEx() const { return ex_; }

  /**
   * @brief Gets number of elements in y direction
   * @return Number of elements along y-axis
   */
  ScalarType GetEy() const { return ey_; }

  /**
   * @brief Gets number of elements in z direction
   * @return Number of elements along z-axis
   */
  ScalarType GetEz() const { return ez_; }

  /**
   * @brief Gets domain length in x direction
   * @return Physical length of domain along x-axis
   */
  FloatType GetLx() const { return lx_; }

  /**
   * @brief Gets domain length in y direction
   * @return Physical length of domain along y-axis
   */
  FloatType GetLy() const { return ly_; }

  /**
   * @brief Gets domain length in z direction
   * @return Physical length of domain along z-axis
   */
  FloatType GetLz() const { return lz_; }

  /**
   * @brief Gets spectral element order
   * @return Order of spectral elements (1-5)
   */
  int GetOrder() const { return order_; }

  /**
   * @brief Checks if density model was loaded
   * @return true if density was provided, false for acoustic models
   */
  bool HasDensity() const { return has_density_; }

 private:
  /**
   * @brief Path to velocity SEP header file
   */
  std::string vp_sep_file_;

  /**
   * @brief Path to density SEP header file (empty for acoustic models)
   */
  std::string rho_sep_file_;

  /**
   * @brief Path to velocity binary data file (extracted from header)
   */
  std::string vp_data_file_path_;

  /**
   * @brief Path to density binary data file (extracted from header)
   */
  std::string rho_data_file_path_;

  /**
   * @brief Spectral element order (1-5)
   */
  int order_;

  /**
   * @brief Flag indicating if model is stored on nodes (true) or elements
   * (false)
   */
  bool is_model_on_nodes_;

  /**
   * @brief Flag indicating if density model is available
   */
  bool has_density_;

  /**
   * @brief Storage for parsed SEP header parameters (key-value pairs)
   */
  std::map<std::string, std::string> header_params_;

  /**
   * @brief Number of elements in x direction
   */
  ScalarType ex_;

  /**
   * @brief Number of elements in y direction
   */
  ScalarType ey_;

  /**
   * @brief Number of elements in z direction
   */
  ScalarType ez_;

  /**
   * @brief Domain length in x direction
   */
  FloatType lx_;

  /**
   * @brief Domain length in y direction
   */
  FloatType ly_;

  /**
   * @brief Domain length in z direction
   */
  FloatType lz_;

  /**
   * @brief Element spacing in x direction
   */
  FloatType dx_;

  /**
   * @brief Element spacing in y direction
   */
  FloatType dy_;

  /**
   * @brief Element spacing in z direction
   */
  FloatType dz_;

  /**
   * @brief Global node indices for each element
   * Shape: [n_element, nodes_per_element]
   */
  ARRAY_INT_VIEW global_node_index_;

  /**
   * @brief X-coordinates of all nodes in mesh
   */
  VECTOR_REAL_VIEW nodes_coords_x_;

  /**
   * @brief Y-coordinates of all nodes in mesh
   */
  VECTOR_REAL_VIEW nodes_coords_y_;

  /**
   * @brief Z-coordinates of all nodes in mesh
   */
  VECTOR_REAL_VIEW nodes_coords_z_;

  /**
   * @brief Velocity model stored on nodes
   */
  VECTOR_REAL_VIEW model_vp_node_;

  /**
   * @brief Velocity model stored on elements
   */
  VECTOR_REAL_VIEW model_vp_element_;

  /**
   * @brief Density model stored on nodes
   */
  VECTOR_REAL_VIEW model_rho_node_;

  /**
   * @brief Density model stored on elements
   */
  VECTOR_REAL_VIEW model_rho_element_;

  /**
   * @brief Parses SEP ASCII header file
   *
   * Reads the SEP header file line by line, extracting key-value pairs.
   * Comments (lines starting with #) are ignored. Quotes around values are
   * removed.
   *
   * @param sep_header_file Path to SEP header file to parse
   *
   * @throws std::runtime_error if file cannot be opened
   *
   * @note Results are stored in header_params_ member variable
   * @note This method clears any existing header parameters before parsing
   */
  void ParseSepHeader(const std::string& sep_header_file) {
    std::ifstream file(sep_header_file);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open SEP header file: " +
                               sep_header_file);
    }

    header_params_.clear();
    std::string line;
    while (std::getline(file, line)) {
      size_t comment_pos = line.find('#');
      if (comment_pos != std::string::npos) {
        line = line.substr(0, comment_pos);
      }

      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);

      if (line.empty()) continue;

      size_t eq_pos = line.find('=');
      if (eq_pos != std::string::npos) {
        std::string key = line.substr(0, eq_pos);
        std::string value = line.substr(eq_pos + 1);

        value.erase(std::remove(value.begin(), value.end(), '\"'),
                    value.end());
        value.erase(std::remove(value.begin(), value.end(), '\''),
                    value.end());

        header_params_[key] = value;
      }
    }
    file.close();
  }

  /**
   * @brief Extracts model parameters from parsed header
   *
   * Reads dimensions (n1, n2, n3), sampling intervals (d1, d2, d3), and data
   * file path from the parsed header parameters. Computes domain lengths from
   * dimensions and sampling intervals.
   *
   * SEP Convention:
   * - n1, n2, n3: number of samples in x, y, z directions
   * - d1, d2, d3: sampling intervals in x, y, z directions
   * - in: path to binary data file
   *
   * @throws std::runtime_error if required parameters n1, n2, d1, or d2 are
   * missing
   *
   * @note If n3 or d3 are missing, defaults to 1 and d2 respectively (2D case)
   * @note If "in" parameter is missing, defaults to header filename + "@"
   * @note Handles relative paths by prepending base directory of header file
   */
  void ExtractModelParameters() {
    if (header_params_.find("n1") != header_params_.end()) {
      ex_ = static_cast<ScalarType>(std::stoull(header_params_["n1"]));
    } else {
      throw std::runtime_error("SEP header missing 'n1' parameter");
    }

    if (header_params_.find("n2") != header_params_.end()) {
      ey_ = static_cast<ScalarType>(std::stoull(header_params_["n2"]));
    } else {
      throw std::runtime_error("SEP header missing 'n2' parameter");
    }

    if (header_params_.find("n3") != header_params_.end()) {
      ez_ = static_cast<ScalarType>(std::stoull(header_params_["n3"]));
    } else {
      ez_ = static_cast<ScalarType>(1);
    }

    if (header_params_.find("d1") != header_params_.end()) {
      dx_ = static_cast<FloatType>(std::stod(header_params_["d1"]));
    } else {
      throw std::runtime_error("SEP header missing 'd1' parameter");
    }

    if (header_params_.find("d2") != header_params_.end()) {
      dy_ = static_cast<FloatType>(std::stod(header_params_["d2"]));
    } else {
      throw std::runtime_error("SEP header missing 'd2' parameter");
    }

    if (header_params_.find("d3") != header_params_.end()) {
      dz_ = static_cast<FloatType>(std::stod(header_params_["d3"]));
    } else {
      dz_ = dy_;
    }

    lx_ = dx_ * ex_;
    ly_ = dy_ * ey_;
    lz_ = dz_ * ez_;

    if (header_params_.find("in") != header_params_.end()) {
      vp_data_file_path_ = header_params_["in"];
      if (vp_data_file_path_[0] != '/') {
        size_t last_slash = vp_sep_file_.find_last_of("/\\");
        if (last_slash != std::string::npos) {
          std::string base_path = vp_sep_file_.substr(0, last_slash + 1);
          vp_data_file_path_ = base_path + vp_data_file_path_;
        }
      }
    } else {
      vp_data_file_path_ = vp_sep_file_ + "@";
    }
  }

  /**
   * @brief Validates that density model geometry matches velocity model
   *
   * Parses the density SEP header and verifies that dimensions (n1, n2, n3)
   * match those of the velocity model. This ensures consistency between
   * physical properties.
   *
   * @throws std::runtime_error if dimensions don't match
   * @throws std::runtime_error if density header file cannot be opened
   *
   * @note Also extracts density data file path from header
   */
  void ValidateDensityGeometry() {
    ParseSepHeader(rho_sep_file_);

    ScalarType rho_ex =
        static_cast<ScalarType>(std::stoull(header_params_.at("n1")));
    ScalarType rho_ey =
        static_cast<ScalarType>(std::stoull(header_params_.at("n2")));
    ScalarType rho_ez = header_params_.find("n3") != header_params_.end()
                            ? static_cast<ScalarType>(
                                  std::stoull(header_params_["n3"]))
                            : static_cast<ScalarType>(1);

    if (rho_ex != ex_ || rho_ey != ey_ || rho_ez != ez_) {
      throw std::runtime_error(
          "Density model dimensions do not match velocity model: vp(" +
          std::to_string(ex_) + "," + std::to_string(ey_) + "," +
          std::to_string(ez_) + ") vs rho(" + std::to_string(rho_ex) + "," +
          std::to_string(rho_ey) + "," + std::to_string(rho_ez) + ")");
    }

    if (header_params_.find("in") != header_params_.end()) {
      rho_data_file_path_ = header_params_["in"];
      if (rho_data_file_path_[0] != '/') {
        size_t last_slash = rho_sep_file_.find_last_of("/\\");
        if (last_slash != std::string::npos) {
          std::string base_path = rho_sep_file_.substr(0, last_slash + 1);
          rho_data_file_path_ = base_path + rho_data_file_path_;
        }
      }
    } else {
      rho_data_file_path_ = rho_sep_file_ + "@";
    }
  }

  /**
   * @brief Initializes global node index array for all elements
   *
   * Creates the connectivity array that maps local node indices within each
   * element to global node indices in the mesh. This accounts for node sharing
   * between adjacent elements in the spectral element method.
   *
   * The indexing follows a structured pattern:
   * - Elements are numbered in i-j-k order
   * - Within each element, nodes are numbered in x-y-z order
   * - Global nodes account for overlapping nodes between elements
   *
   * @note Allocates global_node_index_ with shape [n_element,
   * nodes_per_element]
   * @note nodes_per_element = (order+1)^3
   */
  void InitGlobalNodeList() {
    int nodes_x = order_ + 1;
    int nodes_y = order_ + 1;
    int nodes_z = order_ + 1;
    int total_nodes = nodes_x * nodes_y * nodes_z;
    global_node_index_ = allocateArray2D<ARRAY_INT_VIEW>(
        ex_ * ey_ * ez_, total_nodes, "global node index");
    int nx = ex_ * order_ + 1;
    int ny = ey_ * order_ + 1;
    int nz = ez_ * order_ + 1;

    for (int k = 0; k < ez_; k++) {
      for (int j = 0; j < ey_; j++) {
        for (int i = 0; i < ex_; i++) {
          int element_num = i + j * ex_ + k * ex_ * ey_;
          int offset = i * order_ + j * order_ * nx + k * order_ * nx * ny;

          for (int m = 0; m < order_ + 1; m++) {
            for (int n = 0; n < order_ + 1; n++) {
              for (int l = 0; l < order_ + 1; l++) {
                int dof_local =
                    l + n * (order_ + 1) + m * (order_ + 1) * (order_ + 1);
                int dof_global = offset + l + n * nx + m * nx * ny;
                global_node_index_(element_num, dof_local) = dof_global;
              }
            }
          }
        }
      }
    }
  }

  /**
   * @brief Computes GLL node coordinates in one direction
   *
   * Calculates Gauss-Lobatto-Legendre (GLL) quadrature points for a single
   * element in one spatial direction. These points are optimal for spectral
   * element methods.
   *
   * GLL points used by order:
   * - Order 1: Endpoints only
   * - Order 2: Endpoints + midpoint
   * - Order 3: Endpoints + sqrt(1/5) interior points
   * - Order 4: Endpoints + sqrt(3/7) and 0 interior points
   * - Order 5: Endpoints + 4 interior points based on Legendre roots
   *
   * @param h Element size in this direction
   * @param n_element Element index in this direction
   * @param coord Output array for coordinates (size = order+1)
   *
   * @note GLL points are in reference element [-1, 1], then mapped to physical
   * space
   */
  void GetCoordInOneDirection(const int& h, const int& n_element,
                              float* coord) {
    float xi[MAX_ORDER + 1];

    switch (order_) {
      case 1:
        xi[0] = -1.f;
        xi[1] = 1.f;
        break;
      case 2:
        xi[0] = -1.f;
        xi[1] = 0.f;
        xi[2] = 1.f;
        break;
      case 3: {
        static constexpr float sqrt5 = 2.2360679774997897f;
        xi[0] = -1.0f;
        xi[1] = -1.f / sqrt5;
        xi[2] = 1.f / sqrt5;
        xi[3] = 1.f;
        break;
      }
      case 4: {
        static constexpr float sqrt3_7 = 0.6546536707079771f;
        xi[0] = -1.0f;
        xi[1] = -sqrt3_7;
        xi[2] = 0.0f;
        xi[3] = sqrt3_7;
        xi[4] = 1.0f;
        break;
      }
      case 5: {
        static constexpr float sqrt__7_plus_2sqrt7__ = 3.50592393273573196f;
        static constexpr float sqrt__7_mins_2sqrt7__ = 1.30709501485960033f;
        static constexpr float sqrt_inv21 = 0.218217890235992381f;
        xi[0] = -1.0f;
        xi[1] = -sqrt_inv21 * sqrt__7_plus_2sqrt7__;
        xi[2] = -sqrt_inv21 * sqrt__7_mins_2sqrt7__;
        xi[3] = sqrt_inv21 * sqrt__7_mins_2sqrt7__;
        xi[4] = sqrt_inv21 * sqrt__7_plus_2sqrt7__;
        xi[5] = 1.0f;
        break;
      }
      default:
        break;
    }

    int i = n_element;
    float x0 = i * h;
    float x1 = (i + 1) * h;
    float b = (x1 + x0) / 2.f;
    float a = b - x0;

    for (int j = 0; j < order_ + 1; j++) {
      coord[j] = a * xi[j] + b;
    }
  }

  /**
   * @brief Initializes coordinates for all nodes in the mesh
   *
   * Computes physical coordinates for all GLL nodes in the mesh by:
   * 1. Determining element positions based on mesh dimensions
   * 2. Computing GLL points within each element
   * 3. Mapping to global node indices
   *
   * This creates three separate coordinate arrays (x, y, z) indexed by global
   * node number.
   *
   * @note Allocates nodes_coords_x_, nodes_coords_y_, nodes_coords_z_
   * @note Total nodes = (ex*order+1) * (ey*order+1) * (ez*order+1)
   */
  void InitNodesCoords() {
    int nodes_x = ex_ * order_ + 1;
    int nodes_y = ey_ * order_ + 1;
    int nodes_z = ez_ * order_ + 1;
    int total_nodes = nodes_x * nodes_y * nodes_z;

    nodes_coords_x_ =
        allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords x");
    nodes_coords_y_ =
        allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords y");
    nodes_coords_z_ =
        allocateVector<VECTOR_REAL_VIEW>(total_nodes, "nodes coords z");

    float coord_x[MAX_ORDER + 1];
    float coord_y[MAX_ORDER + 1];
    float coord_z[MAX_ORDER + 1];

    auto hx = lx_ / ex_;
    auto hy = ly_ / ey_;
    auto hz = lz_ / ez_;

    for (int n = 0; n < ez_; n++) {
      GetCoordInOneDirection(hz, n, coord_z);
      for (int m = 0; m < ey_; m++) {
        GetCoordInOneDirection(hy, m, coord_y);
        for (int l = 0; l < ex_; l++) {
          GetCoordInOneDirection(hx, l, coord_x);

          for (int k = 0; k < order_ + 1; k++) {
            for (int j = 0; j < order_ + 1; j++) {
              for (int i = 0; i < order_ + 1; i++) {
                int global_i = l * order_ + i;
                int global_j = m * order_ + j;
                int global_k = n * order_ + k;

                int global_node_index =
                    global_i + global_j * nodes_x + global_k * nodes_x * nodes_y;

                if (global_i < nodes_x && global_j < nodes_y &&
                    global_k < nodes_z) {
                  nodes_coords_x_(global_node_index) = coord_x[i];
                  nodes_coords_y_(global_node_index) = coord_y[j];
                  nodes_coords_z_(global_node_index) = coord_z[k];
                }
              }
            }
          }
        }
      }
    }
  }

  /**
   * @brief Loads velocity model from SEP binary file
   *
   * Reads velocity data from the binary file specified in the SEP header.
   * Supports native_float and native_double formats. Data is read as
   * element-centered values and distributed to either nodes or elements based
   * on is_model_on_nodes_ flag.
   *
   * @throws std::runtime_error if data file cannot be opened
   * @throws std::runtime_error if data format is unsupported
   *
   * @note If is_model_on_nodes_ is true, calls MapElementDataToNodes()
   * @note If is_model_on_nodes_ is false, directly assigns to elements
   */
  void LoadVelocity() {
    int n_element = ex_ * ey_ * ez_;
    int n_node = (ex_ * order_ + 1) * (ey_ * order_ + 1) * (ez_ * order_ + 1);

    vectorReal vp_sep_data;
    size_t total_elements = static_cast<size_t>(ex_) *
                            static_cast<size_t>(ey_) *
                            static_cast<size_t>(ez_);
    vp_sep_data.resize(total_elements);

    std::ifstream file(vp_data_file_path_, std::ios::binary);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open velocity data file: " +
                               vp_data_file_path_);
    }

    size_t element_size = sizeof(float);
    std::string data_format = "native_float";
    if (header_params_.find("esize") != header_params_.end()) {
      element_size = std::stoull(header_params_["esize"]);
    }
    if (header_params_.find("data_format") != header_params_.end()) {
      data_format = header_params_["data_format"];
    }

    if (data_format == "native_float" && element_size == sizeof(float)) {
      std::vector<float> temp(total_elements);
      file.read(reinterpret_cast<char*>(temp.data()),
                total_elements * sizeof(float));
      for (size_t i = 0; i < total_elements; ++i) {
        vp_sep_data[i] = static_cast<Real>(temp[i]);
      }
    } else if (data_format == "native_double" &&
