/**
 * @file validate_solution.cpp
 * @brief Validation system for FUnTiDES spectral element solver
 *
 * This program validates FUnTiDES numerical solutions by comparing them
 * against analytical reference solutions. It supports:
 * - Multiple polynomial orders (1-3)
 * - Different mesh types (cartesian, ucartesian)
 * - Acoustic and elastic physics
 * - Model-on-nodes configuration
 * - Automatic normalization to handle amplitude differences
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @struct ComparisonResult
 * @brief Stores error metrics from solution comparison
 */
struct ComparisonResult
{
  double l2_error;    ///< L2 norm of the error
  double linf_error;  ///< L-infinity norm of the error
  double
      relative_l2_error;  ///< Relative L2 error (normalized by analytical norm)
  size_t num_points;      ///< Number of points compared
};

/**
 * @brief Execute a system command and check its return code
 *
 * @param command The command string to execute
 * @return true if command succeeded (return code 0), false otherwise
 */
bool execute_command(const std::string& command)
{
  std::cout << "Executing: " << command << std::endl;
  int ret = std::system(command.c_str());
  if (ret != 0)
  {
    std::cerr << "Command failed with return code: " << ret << std::endl;
    return false;
  }
  return true;
}

/**
 * @brief Read a specific column from a text file
 *
 * Reads space or tab separated data from a text file, skipping empty lines
 * and comments (lines starting with '#'). Optionally inverts the sign of
 * the values (useful for elastic formulation differences).
 *
 * @param filename Path to the input file
 * @param column_index Zero-based column index to read (0 = first column)
 * @param invert_sign If true, multiply all values by -1
 * @return Vector containing the column values
 * @throw std::runtime_error if file cannot be opened
 */
std::vector<double> read_column(const std::string& filename, int column_index,
                                bool invert_sign = false)
{
  std::vector<double> data;
  std::ifstream file(filename);

  if (!file.is_open())
  {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  std::string line;
  while (std::getline(file, line))
  {
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#')
    {
      continue;
    }

    std::istringstream iss(line);
    std::vector<double> row_values;
    double value;

    while (iss >> value)
    {
      row_values.push_back(value);
    }

    if (column_index < static_cast<int>(row_values.size()))
    {
      double column_value = row_values[column_index];
      // Invert sign if requested (for elastic case where formulations differ)
      if (invert_sign)
      {
        column_value = -column_value;
      }
      data.push_back(column_value);
    }
  }

  return data;
}

/**
 * @brief Compare numerical and analytical solutions and compute error metrics
 *
 * Compares two solution vectors point-by-point and computes various error
 * metrics. Optionally normalizes both solutions to unit norm before comparison,
 * which makes the comparison amplitude-independent and focuses on the solution
 * shape rather than absolute values.
 *
 * @param numerical Numerical solution from FUnTiDES
 * @param analytical Analytical reference solution
 * @param normalize If true, normalize both solutions before comparison
 * (default: true)
 * @return ComparisonResult structure containing error metrics
 * @throw std::runtime_error if solution sizes don't match
 */
ComparisonResult compare_solutions(const std::vector<double>& numerical,
                                   const std::vector<double>& analytical,
                                   bool normalize = true)
{
  if (numerical.size() != analytical.size())
  {
    throw std::runtime_error(
        "Size mismatch: numerical (" + std::to_string(numerical.size()) +
        ") vs analytical (" + std::to_string(analytical.size()) + ")");
  }

  // Make copies for normalization
  std::vector<double> num_data = numerical;
  std::vector<double> ana_data = analytical;

  // Normalize if requested (to handle amplitude differences)
  if (normalize)
  {
    // Calculate norms
    double num_norm = 0.0;
    double ana_norm = 0.0;
    for (size_t i = 0; i < num_data.size(); ++i)
    {
      num_norm += num_data[i] * num_data[i];
      ana_norm += ana_data[i] * ana_data[i];
    }
    num_norm = std::sqrt(num_norm);
    ana_norm = std::sqrt(ana_norm);

    // Normalize both solutions to unit norm
    if (num_norm > 1e-14 && ana_norm > 1e-14)
    {
      for (size_t i = 0; i < num_data.size(); ++i)
      {
        num_data[i] /= num_norm;
        ana_data[i] /= ana_norm;
      }
      std::cout << "\nNormalization applied:" << std::endl;
      std::cout << "  Original ||numerical||_2 = " << num_norm << std::endl;
      std::cout << "  Original ||analytical||_2 = " << ana_norm << std::endl;
      std::cout << "  Amplitude ratio: " << num_norm / ana_norm << std::endl;
    }
  }

  ComparisonResult result;
  result.num_points = num_data.size();
  result.l2_error = 0.0;
  result.linf_error = 0.0;
  double analytical_norm = 0.0;

  size_t max_error_index = 0;

  // Compute error metrics
  for (size_t i = 0; i < num_data.size(); ++i)
  {
    double diff = std::abs(num_data[i] - ana_data[i]);
    result.l2_error += diff * diff;
    if (diff > result.linf_error)
    {
      result.linf_error = diff;
      max_error_index = i;
    }
    analytical_norm += ana_data[i] * ana_data[i];
  }

  result.l2_error = std::sqrt(result.l2_error);
  analytical_norm = std::sqrt(analytical_norm);

  // Compute relative error (avoid division by zero)
  if (analytical_norm > 1e-14)
  {
    result.relative_l2_error = result.l2_error / analytical_norm;
  }
  else
  {
    result.relative_l2_error =
        result.l2_error;  // Absolute error if norm is too small
  }

  // Print location of maximum error for debugging
  std::cout << "\nMaximum error location:" << std::endl;
  std::cout << "  Point index:           " << max_error_index << std::endl;
  std::cout << "  Numerical value:       " << num_data[max_error_index]
            << std::endl;
  std::cout << "  Analytical value:      " << ana_data[max_error_index]
            << std::endl;
  std::cout << "  Absolute difference:   " << result.linf_error << std::endl;

  return result;
}

/**
 * @brief Print validation results in a formatted table
 *
 * @param result The comparison result to display
 * @param tolerance The tolerance threshold for validation
 */
void print_results(const ComparisonResult& result, double tolerance)
{
  std::cout << "\n" << std::string(60, '=') << std::endl;
  std::cout << "VALIDATION RESULTS" << std::endl;
  std::cout << std::string(60, '=') << std::endl;
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "Number of points:     " << result.num_points << std::endl;
  std::cout << "L2 error:             " << result.l2_error << std::endl;
  std::cout << "L_inf error:          " << result.linf_error << std::endl;
  std::cout << "Relative L2 error:    " << result.relative_l2_error
            << std::endl;
  std::cout << "Tolerance:            " << tolerance << std::endl;
  std::cout << std::string(60, '=') << std::endl;

  if (result.relative_l2_error < tolerance)
  {
    std::cout << "✓ VALIDATION PASSED" << std::endl;
  }
  else
  {
    std::cout << "✗ VALIDATION FAILED" << std::endl;
  }
  std::cout << std::string(60, '=') << std::endl;
}

/**
 * @brief Main validation program
 *
 * Workflow:
 * 1. Parse command line arguments
 * 2. Run FUnTiDES solver with specified parameters
 * 3. Post-process results to extract receiver data
 * 4. Compare with analytical reference solution
 * 5. Report validation results
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 * @return 0 if validation passed, 1 if validation failed or error occurred
 */
int main(int argc, char* argv[])
{
  // ====================================================================
  // Default paths (can be overridden by CMake at compile time)
  // ====================================================================

#ifndef FUNTIDES_EXECUTABLE
// Fallback if not defined by CMake (for manual compilation)
#define FUNTIDES_EXECUTABLE "./bin/funtides-sem"
#endif

#ifndef EXTRACT_RECEIVERS_SCRIPT
#define EXTRACT_RECEIVERS_SCRIPT "../scripts/adios/adios_single_receiver_viz.py"
#endif

#ifndef VALIDATION_REFERENCE_ACOUSTIC
#define VALIDATION_REFERENCE_ACOUSTIC \
  "../tests/validation/analyticalsolution/P.dat"
#endif

#ifndef VALIDATION_REFERENCE_ELASTIC
#define VALIDATION_REFERENCE_ELASTIC \
  "../tests/validation/analyticalsolution/Ux.dat"
#endif

  // ====================================================================
  // Configuration variables
  // ====================================================================

  // solver_command will be built after parsing command line arguments
  std::string solver_command;

  // Set MPLBACKEND=Agg to prevent matplotlib from showing plots
  // (non-interactive mode)
  std::string postprocess_command =
      std::string("MPLBACKEND=Agg python3 ") + EXTRACT_RECEIVERS_SCRIPT;

  // Files to compare
  std::string numerical_file = "seismogram_output.txt";  // FUnTiDES output
  std::string analytical_file;  // Will be set based on test type

  // Column indices for comparison (0-based indexing)
  int numerical_column = 0;   // For FUnTiDES output, solution is in 1st column
  int analytical_column = 1;  // For reference files, solution is in 2nd column

  // Validation tolerance (relative L2 error threshold)
  double tolerance = 5e-2;  // Default 5%, usually overridden by CMake

  // ====================================================================
  // Parse command line arguments (optional override)
  // ====================================================================

  int polynomial_order = 1;             // Default polynomial order
  bool is_elastic = false;              // Default: acoustic case
  bool is_model_on_nodes = false;       // Default: model not on nodes
  std::string mesh_type = "cartesian";  // Default: cartesian mesh

  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (arg == "--tolerance" && i + 1 < argc)
    {
      tolerance = std::stod(argv[++i]);
    }
    else if (arg == "--order" && i + 1 < argc)
    {
      polynomial_order = std::stoi(argv[++i]);
    }
    else if (arg == "--elastic" || arg == "--is-elastic")
    {
      is_elastic = true;
    }
    else if (arg == "--is-model-on-nodes")
    {
      is_model_on_nodes = true;
    }
    else if (arg == "--mesh" && i + 1 < argc)
    {
      mesh_type = argv[++i];
    }
    else if (arg == "--help")
    {
      std::cout << "Usage: " << argv[0] << " [OPTIONS]" << std::endl;
      std::cout << "Options:" << std::endl;
      std::cout << "  --tolerance VALUE    Validation tolerance (default: "
                << tolerance << ")" << std::endl;
      std::cout << "  --order N            Polynomial order (default: 1)"
                << std::endl;
      std::cout << "  --elastic            Run elastic case (default: acoustic)"
                << std::endl;
      std::cout << "  --is-model-on-nodes  Enable model on nodes (default: off)"
                << std::endl;
      std::cout << "  --mesh TYPE          Mesh type: cartesian or ucartesian "
                   "(default: cartesian)"
                << std::endl;
      std::cout << "  --help               Show this help message" << std::endl;
      return 0;
    }
  }

  // ====================================================================
  // Build solver command with parameters
  // ====================================================================

  // Build optional flags
  std::string elastic_flag = is_elastic ? " --is-elastic" : "";
  std::string model_on_nodes_flag =
      is_model_on_nodes ? " --is-model-on-nodes" : "";

  // Adjust mesh resolution based on polynomial order to keep similar number of
  // DOFs For order p, an element has (p+1) points per direction To keep total
  // points constant: n_elem * (p+1) = constant Base: order 1 with 200 elements
  // → 200*2 = 400 points per direction
  int base_points = 400;  // Target number of points per direction
  int elements_per_direction = base_points / (polynomial_order + 1);

  // Construct full solver command
  solver_command = std::string(FUNTIDES_EXECUTABLE) + " --order " +
                   std::to_string(polynomial_order) + " --mesh=" + mesh_type +
                   " --ex " + std::to_string(elements_per_direction) +
                   " --ey " + std::to_string(elements_per_direction) +
                   " --ez " + std::to_string(elements_per_direction) +
                   " --implem makutu --dt 0.001 --anisotropy iso " +
                   elastic_flag + model_on_nodes_flag;

  // Select reference file based on test type
  std::string test_type = is_elastic ? "elastic" : "acoustic";
  analytical_file =
      is_elastic ? VALIDATION_REFERENCE_ELASTIC : VALIDATION_REFERENCE_ACOUSTIC;

  // ====================================================================
  // Display configuration
  // ====================================================================

  std::cout << "\n" << std::string(60, '=') << std::endl;
  std::cout << "VALIDATION CONFIGURATION" << std::endl;
  std::cout << std::string(60, '=') << std::endl;
  std::cout << "Test type:            " << test_type << std::endl;
  std::cout << "Polynomial order:     " << polynomial_order << std::endl;
  std::cout << "Mesh type:            " << mesh_type << std::endl;
  std::cout << "Model on nodes:       " << (is_model_on_nodes ? "yes" : "no")
            << std::endl;
  std::cout << "Mesh resolution:      " << elements_per_direction << "x"
            << elements_per_direction << "x" << elements_per_direction
            << " elements" << std::endl;
  std::cout << "Points per direction: ~"
            << elements_per_direction * (polynomial_order + 1) << std::endl;
  std::cout << "Tolerance:            " << tolerance << std::endl;
  std::cout << "Reference file:       " << analytical_file << std::endl;
  std::cout << std::string(60, '=') << std::endl;

  // ====================================================================
  // EXECUTION
  // ====================================================================

  try
  {
    // Step 1: Run the solver
    std::cout << "\n[1/3] Running FUnTiDES solver..." << std::endl;
    if (!execute_command(solver_command))
    {
      std::cerr << "Solver execution failed!" << std::endl;
      return 1;
    }

    // Step 2: Post-process results to generate receiver file
    std::cout << "\n[2/3] Post-processing results..." << std::endl;
    if (!execute_command(postprocess_command))
    {
      std::cerr << "Post-processing failed!" << std::endl;
      return 1;
    }

    // Step 3: Compare results
    std::cout << "\n[3/3] Comparing solutions..." << std::endl;

    // Read numerical solution
    std::vector<double> numerical =
        read_column(numerical_file, numerical_column);

    // Read analytical solution (invert sign for elastic case due to formulation
    // differences)
    std::vector<double> analytical =
        read_column(analytical_file, analytical_column, is_elastic);

    std::cout << "Read " << numerical.size()
              << " points from numerical solution" << std::endl;
    std::cout << "Read " << analytical.size()
              << " points from reference solution";
    if (is_elastic)
    {
      std::cout << " (sign inverted for elastic formulation)";
    }
    std::cout << std::endl;

    // Compare solutions (with automatic normalization)
    ComparisonResult result = compare_solutions(numerical, analytical);

    // Debug: show first few values to check alignment
    std::cout << "\nFirst 5 values comparison:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), numerical.size()); ++i)
    {
      std::cout << "  Point " << i << ": numerical=" << numerical[i]
                << ", analytical=" << analytical[i]
                << ", diff=" << std::abs(numerical[i] - analytical[i])
                << std::endl;
    }

    // Print results and determine pass/fail
    print_results(result, tolerance);

    // Return appropriate exit code
    return (result.relative_l2_error < tolerance) ? 0 : 1;
  }
  catch (const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
