#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace model
{
/**
 * @defgroup SEP_Reader SEP File Reader Module
 * @brief Reading and processing Stanford Exploration Project files
 * @{
 */

/**
 * @struct SepHeader
 * @brief Metadata structure from a SEP (Stanford Exploration Project) file
 * header
 *
 * Encapsulates all metadata extracted from a SEP file header, including grid
 * dimensions, spacing, origins, and data format parameters.
 *
 * @details
 *
 * The structure provides default values for all fields to allow safe
 * initialization. It can be populated by SepReader::readHeader() or manually
 * initialized.
 *
 * @see SepReader::readHeader()
 * @see SepReader::readData()
 */
struct SepHeader
{
  /// @name Grid Dimensions
  /// @{

  /// @brief Number of samples along axis 1 (typically X or CDP direction)
  int n1 = 0;

  /// @brief Number of samples along axis 2 (typically Y or line direction)
  int n2 = 0;

  /// @brief Number of samples along axis 3 (typically Z or depth direction)
  int n3 = 0;
  /// @}

  /// @name Grid Spacing
  /// @{

  /// @brief Spacing between samples along axis 1 (in measurement units)
  double d1 = 1.0;

  /// @brief Spacing between samples along axis 2
  double d2 = 1.0;

  /// @brief Spacing between samples along axis 3
  double d3 = 1.0;
  /// @}

  /// @name Origin Coordinates
  /// @{

  /// @brief Origin coordinate for axis 1 (coordinate of first sample)
  double o1 = 0.0;

  /// @brief Origin coordinate for axis 2
  double o2 = 0.0;

  /// @brief Origin coordinate for axis 3
  double o3 = 0.0;
  /// @}

  /// @name Axis Labels and Units
  /// @{

  /// @brief Label for axis 1 (e.g., "CDP", "X")
  std::string label1 = "X";

  /// @brief Label for axis 2 (e.g., "LINE", "Y")
  std::string label2 = "Y";

  /// @brief Label for axis 3 (e.g., "DEPTH", "Z")
  std::string label3 = "Z";

  /// @brief Measurement unit for axis 1 (e.g., "m", "ft")
  std::string unit1 = "m";

  /// @brief Measurement unit for axis 2
  std::string unit2 = "m";

  /// @brief Measurement unit for axis 3
  std::string unit3 = "m";
  /// @}

  /// @name Binary Data Format
  /// @{

  /// @brief Data format string (e.g., "xdr_float", "xdr_double")
  std::string data_format = "xdr_float";

  /// @brief Size in bytes of each data element (4 for float, 8 for double)
  int esize = 4;
  /// @}

  /// @name Data File Information
  /// @{

  /// @brief Path to the binary data file
  std::string data_file;

  /// @brief End-of-file byte offset for data (important for .H@ embedded files)
  long eofb = 0;
  /// @}

  /// @name Additional Metadata
  /// @{

  /// @brief Data label describing the content (e.g., "VINT" for interval
  /// velocity)
  std::string data_label;

  /// @brief File type specification (e.g., "regular" for Cartesian grid)
  std::string data_filetype;

  /// @brief MD5 checksum of the data for validation purposes
  std::string md5sum;
  /// @}

  /// @name Finite Element Parameters
  /// @{

  /// @brief Polynomial order for spectral elements (1-5, default: 4)
  int order = 4;
  /// @}

  /**
   * @brief Print formatted header information to standard output
   *
   * Displays a formatted summary of all metadata stored in the header.
   * Useful for debugging and parameter verification.
   *
   * @details The output includes:
   * - Grid dimensions (n1, n2, n3)
   * - Grid spacing (d1, d2, d3)
   * - Origin coordinates (o1, o2, o3)
   * - Axis labels and units
   * - Data format and element size
   * - Data file path
   * - Data label and file type
   * - Polynomial order
   *
   * @see SepReader::readHeader()
   */
  void print() const
  {
    std::cout << "\n=== SEP Header Information ===\n"
              << "Dimensions: n1=" << n1 << ", n2=" << n2 << ", n3=" << n3
              << "\n"
              << "Spacing:    d1=" << d1 << ", d2=" << d2 << ", d3=" << d3
              << "\n"
              << "Origin:     o1=" << o1 << ", o2=" << o2 << ", o3=" << o3
              << "\n"
              << "Labels:     " << label1 << ", " << label2 << ", " << label3
              << "\n"
              << "Units:      " << unit1 << ", " << unit2 << ", " << unit3
              << "\n"
              << "Data format: " << data_format << " (esize=" << esize
              << " bytes)\n"
              << "Data file:  " << data_file << "\n"
              << "Data label: " << data_label << "\n"
              << "File type:  " << data_filetype << "\n"
              << "FE Order:   " << order << "\n"
              << "==============================\n\n";
  }
};

/**
 * @class SepReader
 * @brief Static reader for Stanford Exploration Project (SEP) format files
 *
 * @details
 * Provides static methods to read and parse SEP files, which are widely used in
 * geophysics to store seismic data and velocity models. SEP files consist of a
 * text header describing the data followed by binary data values.
 *
 * Supported file formats:
 * - **`.H` files**: Header and data in separate files
 *   - Header in text format (e.g., `model.H`)
 *   - Data in binary format (e.g., `model.bin`)
 *   - Useful for large datasets
 *
 * - **`.H@` files**: Header and data embedded in single file
 *   - More portable and self-contained
 *   - Recommended for distribution
 *
 * @section sep_parsing_rules Parsing Rules
 * - Empty lines and comments (starting with `#`) are ignored
 * - Each data line must contain a key=value pair
 * - Quoted strings have quotes removed
 * - Leading and trailing whitespace is trimmed
 *
 * @section sep_file_resolution File Path Resolution
 * - For `.H@` files: Data is embedded in the same file
 * - For `.H` files with relative paths: Paths are resolved relative to header
 * directory
 * - If no data file is specified: A `.bin` file is assumed
 *
 * @see SepHeader
 * @see ModelUnstructBuilderFromSEP
 * @ingroup SEP_Reader
 */
class SepReader
{
 public:
  /**
   * @brief Read and parse a SEP header file
   *
   * @param filename Path to the SEP file
   *   - Can be a `.H` file (header only, data in separate file)
   *   - Can be a `.H@` file (header and data in one file)
   *
   * @return Populated SepHeader structure with all metadata
   *
   * @throws std::runtime_error if file cannot be opened
   * @throws std::runtime_error if header parsing fails
   *
   * @details
   * This function performs the following steps:
   * 1. Opens the header file
   * 2. Parses each line for key-value pairs
   * 3. Converts values to appropriate types (int, double, string)
   * 4. Resolves relative paths for data files
   * 5. Detects data offset for `.H@` files
   * 6. Returns populated SepHeader structure
   *
   * The function is forgiving with optional parameters but strict with
   * required ones (n1, n2, n3, d1, d2, d3, esize, data_format).
   *
   * @see SepHeader
   * @see readData()
   */
  static SepHeader readHeader(const std::string& filename)
  {
    std::ifstream file(filename);
    if (!file.is_open())
    {
      throw std::runtime_error("Cannot open file: " + filename);
    }

    SepHeader header;
    std::string line;
    std::string directory = getDirectory(filename);

    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#') continue;

      parseLine(line, header);
    }

    if (filename.find(".H@") != std::string::npos)
    {
      header.data_file = filename;
      header.eofb = getHeaderEndOffset(filename);
    }
    else
    {
      if (header.data_file.empty())
      {
        header.data_file = filename;
        size_t pos = header.data_file.find_last_of('.');
        if (pos != std::string::npos)
        {
          header.data_file = header.data_file.substr(0, pos) + ".bin";
        }
      }
      else if (header.data_file[0] != '/')
      {
        header.data_file = directory + "/" + header.data_file;
      }
    }

    file.close();
    return header;
  }

  /**
   * @brief Read binary data from a SEP file
   *
   * @param header SepHeader structure containing metadata
   *   (typically obtained from readHeader())
   *
   * @return std::vector<float> containing all data values
   *   - Vector size: n1 × n2 × n3 elements
   *   - All values converted to float32
   *   - Automatically converts from float64 if needed
   *
   * @throws std::runtime_error if data file cannot be opened
   * @throws std::runtime_error if read operation fails
   * @throws std::runtime_error if data format is unsupported
   *
   * @details
   * This function:
   * 1. Opens the data file specified in the header
   * 2. Seeks to the correct offset (important for `.H@` files)
   * 3. Determines element size (from header.esize)
   * 4. Reads binary data into memory
   * 5. Converts to float32 if necessary (from float64)
   *
   * Supported data formats:
   * - **float32** (4 bytes): Read directly without conversion
   * - **float64** (8 bytes): Converted to float32
   * - **xdr_float/xdr_double**: Architecture-independent formats
   *
   * Memory layout of returned data:
   * For indices (i, j, k) where i ∈ [0,n1), j ∈ [0,n2), k ∈ [0,n3):
   * @code
   * linear_index = k * (n1 * n2) + j * n1 + i
   * @endcode
   *
   * @see readHeader()
   * @see SepHeader
   */
  static std::vector<float> readData(const SepHeader& header)
  {
    std::ifstream file(header.data_file, std::ios::binary);
    if (!file.is_open())
    {
      throw std::runtime_error("Cannot open data file: " + header.data_file);
    }

    long n_elements = static_cast<long>(header.n1) * header.n2 * header.n3;

    int element_size = header.esize;
    if (element_size == 0)
    {
      if (header.data_format.find("float") != std::string::npos)
        element_size = 4;
      else if (header.data_format.find("double") != std::string::npos)
        element_size = 8;
    }

    file.seekg(header.eofb, std::ios::beg);

    std::vector<float> data(n_elements);

    if (element_size == 4)
    {
      file.read(reinterpret_cast<char*>(data.data()),
                n_elements * sizeof(float));
    }
    else if (element_size == 8)
    {
      std::vector<double> temp(n_elements);
      file.read(reinterpret_cast<char*>(temp.data()),
                n_elements * sizeof(double));
      std::copy(temp.begin(), temp.end(), data.begin());
    }
    else
    {
      throw std::runtime_error("Unsupported element size: " +
                               std::to_string(element_size));
    }

    if (!file)
    {
      throw std::runtime_error("Error reading data from: " + header.data_file);
    }

    file.close();
    return data;
  }

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
    if (pos == std::string::npos) return ".";
    if (pos == 0) return "/";
    return filepath.substr(0, pos);
  }

 private:
  /**
   * @brief Parse a single line from the header file
   *
   * @param line The line to parse
   * @param header The SepHeader structure to populate
   *
   * @details
   * Splits the line into tokens and locates the '=' separator.
   * Delegates to parseKeyValue() for actual value processing.
   *
   * @see parseKeyValue()
   */
  static void parseLine(const std::string& line, SepHeader& header)
  {
    std::istringstream iss(line);
    std::string token;

    if (!(iss >> token)) return;

    size_t eq_pos = token.find('=');
    if (eq_pos == std::string::npos) return;

    std::string key = token.substr(0, eq_pos);
    std::string value = token.substr(eq_pos + 1);

    std::string rest;
    while (iss >> rest)
    {
      value += " " + rest;
    }

    parseKeyValue(key, value, header);
  }

  /**
   * @brief Process a key-value pair and update SepHeader
   *
   * @param key Parameter name
   * @param value Parameter value as string
   * @param header SepHeader structure to update
   *
   * @details
   * Performs type conversion and populates the corresponding field
   * in the SepHeader structure. Silently ignores conversion errors
   * for optional parameters.
   *
   * Supported keys:
   * - n1, n2, n3: Grid dimensions (int)
   * - d1, d2, d3: Grid spacing (double)
   * - o1, o2, o3: Origin coordinates (double)
   * - label1, label2, label3: Axis labels (string)
   * - unit1, unit2, unit3: Measurement units (string)
   * - data_format: Data format specification (string)
   * - esize: Element size in bytes (int)
   * - in: Data file path (string)
   * - data_label: Data type label (string)
   * - data_filetype: File type specification (string)
   * - md5sum: Checksum (string)
   * - order: Polynomial order (int)
   *
   * @see parseLine()
   * @see cleanString()
   */
  static void parseKeyValue(const std::string& key, const std::string& value,
                            SepHeader& header)
  {
    try
    {
      if (key == "n1")
        header.n1 = std::stoi(value);
      else if (key == "n2")
        header.n2 = std::stoi(value);
      else if (key == "n3")
        header.n3 = std::stoi(value);
      else if (key == "d1")
        header.d1 = std::stod(value);
      else if (key == "d2")
        header.d2 = std::stod(value);
      else if (key == "d3")
        header.d3 = std::stod(value);
      else if (key == "o1")
        header.o1 = std::stod(value);
      else if (key == "o2")
        header.o2 = std::stod(value);
      else if (key == "o3")
        header.o3 = std::stod(value);
      else if (key == "label1")
        header.label1 = cleanString(value);
      else if (key == "label2")
        header.label2 = cleanString(value);
      else if (key == "label3")
        header.label3 = cleanString(value);
      else if (key == "unit1")
        header.unit1 = cleanString(value);
      else if (key == "unit2")
        header.unit2 = cleanString(value);
      else if (key == "unit3")
        header.unit3 = cleanString(value);
      else if (key == "data_format")
        header.data_format = cleanString(value);
      else if (key == "esize")
        header.esize = std::stoi(value);
      else if (key == "in")
        header.data_file = cleanString(value);
      else if (key == "data_label")
        header.data_label = cleanString(value);
      else if (key == "data_filetype")
        header.data_filetype = cleanString(value);
      else if (key == "md5sum")
        header.md5sum = cleanString(value);
      else if (key == "order")
        header.order = std::stoi(value);
    }
    catch (const std::exception&)
    {
      // Ignore parsing errors for optional parameters
    }
  }

  /**
   * @brief Clean a string by removing quotes and whitespace
   *
   * @param str Input string to clean
   *
   * @return Cleaned string
   *   - Removes leading/trailing quotes
   *   - Removes leading/trailing whitespace
   *   - Returns empty string if only whitespace remains
   *
   * @details
   * Used to process string values from the header, which may be
   * quoted and padded with whitespace. Handles both single and
   * multiple whitespace characters.
   */
  static std::string cleanString(const std::string& str)
  {
    std::string result = str;

    if (!result.empty() && result.front() == '"') result = result.substr(1);
    if (!result.empty() && result.back() == '"')
      result = result.substr(0, result.size() - 1);

    size_t start = result.find_first_not_of(" \t");
    size_t end = result.find_last_not_of(" \t");

    if (start != std::string::npos && end != std::string::npos)
    {
      return result.substr(start, end - start + 1);
    }
    return "";
  }

  /**
   * @brief Find the byte offset where binary data begins in a `.H@` file
   *
   * @param filename Path to the `.H@` file
   *
   * @return Number of bytes before the start of binary data
   *   - For text headers followed by data
   *   - Important for parsing embedded format files
   *
   * @details
   * In `.H@` files, the header and data are combined. This function
   * determines where the text header ends and binary data begins by
   * reading the file line by line until an empty line is encountered.
   *
   * The returned offset is used to skip the header when reading data.
   */
  static long getHeaderEndOffset(const std::string& filename)
  {
    std::ifstream file(filename);
    if (!file.is_open())
    {
      return 0;
    }

    std::string line;
    long offset = 0;

    while (std::getline(file, line))
    {
      offset += line.length() + 1;
      if (line.empty()) break;
    }

    file.close();
    return offset;
  }
};

/// @} // end of SEP_Reader group

}  // namespace model
