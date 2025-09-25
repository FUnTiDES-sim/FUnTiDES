#ifndef ADIOS_UTILS_H_
#define ADIOS_UTILS_H_

#include <adios2.h>

#include <string>

#include "data_type.h"

namespace io
{

template <typename Data>
void write1dArray(adios2::ADIOS& adios, Data& array,
                  const std::string& filename, const std::string dataname)
{
  adios2::IO io = adios.DeclareIO(dataname);

  using ValueType = typename Data::value_type;

  std::vector<size_t> dims = {array.extent(0)};  // For 1D array
  auto var_data = io.DefineVariable<ValueType>(dataname, {}, {}, dims);

  adios2::Engine writer = io.Open(filename, adios2::Mode::Write);

  writer.BeginStep();
  writer.Put(var_data, array.data());
  writer.EndStep();
  writer.Close();
}

template <typename Data>
void write2dArray(adios2::ADIOS& adios, Data array, const std::string& filename,
                  const std::string dataname)
{
  adios2::IO io = adios.DeclareIO(dataname);

  using ValueType = typename Data::value_type;

  std::vector<size_t> dims = {array.extent(0),
                              array.extent(1)};  // For 1D array
  auto var_data = io.DefineVariable<ValueType>(dataname, {}, {}, dims);

  adios2::Engine writer = io.Open(filename, adios2::Mode::Write);

  writer.BeginStep();
  writer.Put(var_data, array.data());
  writer.EndStep();
  writer.Close();
}

}  // namespace io

#endif  // ADIOS_UTILS_H_
