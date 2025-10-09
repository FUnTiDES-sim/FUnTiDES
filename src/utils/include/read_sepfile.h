#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

class SEPReader
{
 public:
  struct Header
  {
    int n1, n2, n3;        // dimensions
    float d1, d2, d3;      // spacing
    float o1, o2, o3;      // origins
    std::string datatype;  // data type
    std::string in;        // binary file path
  };

  static Header readHeader(const std::string& headerFile)
  {
    Header hdr = {};
    std::ifstream file(headerFile);
    if (!file)
    {
      throw std::runtime_error("Cannot open SEP header file: " + headerFile);
    }

    std::string line;
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#') continue;

      std::string key, value;
      size_t pos = line.find('=');
      if (pos != std::string::npos)
      {
        key = line.substr(0, pos);
        value = line.substr(pos + 1);

        // Trim whitespace
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);

        if (key == "n1")
          hdr.n1 = std::stoi(value);
        else if (key == "n2")
          hdr.n2 = std::stoi(value);
        else if (key == "n3")
          hdr.n3 = std::stoi(value);
        else if (key == "d1")
          hdr.d1 = std::stof(value);
        else if (key == "d2")
          hdr.d2 = std::stof(value);
        else if (key == "d3")
          hdr.d3 = std::stof(value);
        else if (key == "o1")
          hdr.o1 = std::stof(value);
        else if (key == "o2")
          hdr.o2 = std::stof(value);
        else if (key == "o3")
          hdr.o3 = std::stof(value);
        else if (key == "data_format")
          hdr.datatype = value;
        else if (key == "in")
          hdr.in = value;
      }
    }
    return hdr;
  }

  static std::vector<float> readBinary(const Header& hdr)
  {
    std::ifstream file(hdr.in, std::ios::binary);
    if (!file)
    {
      throw std::runtime_error("Cannot open SEP data file: " + hdr.in);
    }

    size_t nSamples = hdr.n1 * hdr.n2 * hdr.n3;
    std::vector<float> data(nSamples);

    file.read(reinterpret_cast<char*>(data.data()), nSamples * sizeof(float));

    if (file.gcount() != nSamples * sizeof(float))
    {
      throw std::runtime_error("Failed to read complete SEP data file");
    }

    return data;
  }
};