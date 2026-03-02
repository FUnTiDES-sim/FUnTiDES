#ifndef FUNTIDES_UTILS_INCLUDE_DATA_TYPE_VECTOR_H_
#define FUNTIDES_UTILS_INCLUDE_DATA_TYPE_VECTOR_H_
#include <memory>
#include <vector>

template <class T>
class Vector
{
 public:
  Vector(int numRows) : data_(std::make_shared<std::vector<T>>(numRows)) {}
  Vector() : data_(std::make_shared<std::vector<T>>(0)) {}

  // Iterator range constructor for copying from other containers
  template <typename Iterator>
  Vector(Iterator first, Iterator last)
      : data_(std::make_shared<std::vector<T>>(first, last))
  {
  }

  // Element access
  T &operator()(int index) { return (*data_)[index]; }
  const T &operator()(int index) const { return (*data_)[index]; }
  T &operator[](int index) { return (*data_)[index]; }
  const T &operator[](int index) const { return (*data_)[index]; }

  // Assignment operator (shallow copy)
  Vector &operator=(const Vector &other)
  {
    if (this != &other)
    {
      data_ = other.data_;
    }
    return *this;
  }

  // Size/extent
  size_t extent(int dim) const { return data_->size(); }

  // Data access
  T *data() { return data_->data(); }
  const T *data() const { return data_->data(); }

 private:
  std::shared_ptr<std::vector<T>> data_;
};

template <class T>
class Array2D
{
 public:
  Array2D(int numRows, int numCols)
      : data(std::make_shared<std::vector<std::vector<T>>>(
            numRows, std::vector<T>(numCols, 0)))
  {
  }
  Array2D()
      : data(
            std::make_shared<std::vector<std::vector<T>>>(0, std::vector<T>(0)))
  {
  }

  std::vector<T> &operator[](int index) { return (*data)[index]; }
  const std::vector<T> &operator[](int index) const { return (*data)[index]; }
  T &operator()(int row, int col) { return (*data)[row][col]; }
  T &operator()(int row, int col) const
  {
    return const_cast<T &>((*data)[row][col]);
  }

  // Assignment operator (shallow copy)
  Array2D &operator=(const Array2D &other)
  {
    if (this != &other)
    {
      data = other.data;
    }
    return *this;
  }

  size_t extent(int dim) const
  {
    if (dim == 0) return data->size();
    if (dim == 1 && !data->empty()) return (*data)[0].size();
    return 0;
  }

  std::vector<T> getColumn(int colIndex) const
  {
    if (data->empty() || colIndex >= (*data)[0].size())
    {
      return {};  // Empty vector
    }

    std::vector<T> column(data->size());
    for (size_t i = 0; i < data->size(); ++i)
    {
      column[i] = (*data)[i][colIndex];
    }
    return column;
  }

 private:
  std::shared_ptr<std::vector<std::vector<T>>> data;
};

template <class T>
class Array3D
{
 public:
  Array3D(int X, int Y, int Z)
      : data(std::make_shared<std::vector<std::vector<std::vector<T>>>>(
            X, std::vector<std::vector<T>>(Y, std::vector<T>(Z))))
  {
  }
  Array3D()
      : data(std::make_shared<std::vector<std::vector<std::vector<T>>>>(
            0, std::vector<std::vector<T>>(0)))
  {
  }

  std::vector<std::vector<T>> &operator[](int index) { return (*data)[index]; }
  T &operator()(size_t X, size_t Y, size_t Z) { return (*data)[X][Y][Z]; }
  const T &operator()(size_t X, size_t Y, size_t Z) const
  {
    return (*data)[X][Y][Z];
  }

  // Assignment operator (shallow copy)
  Array3D &operator=(const Array3D &other)
  {
    if (this != &other)
    {
      data = other.data;
    }
    return *this;
  }

  size_t extent(int dim) const
  {
    if (dim == 0) return data->size();
    if (dim == 1 && !data->empty()) return (*data)[0].size();
    if (dim == 2 && !data->empty() && !(*data)[0].empty())
      return (*data)[0][0].size();
    return 0;
  }

 private:
  std::shared_ptr<std::vector<std::vector<std::vector<T>>>> data;
};

using vectorReal = Vector<float>;
using vectorInt = Vector<int>;
using vectorDouble = Vector<double>;

using arrayInt = Array2D<int>;
using arrayReal = Array2D<float>;
using arrayDouble = Array2D<double>;

using array3DInt = Array3D<int>;
using array3DReal = Array3D<float>;
using array3DDouble = Array3D<double>;
#endif  // FUNTIDES_UTILS_INCLUDE_DATA_TYPE_VECTOR_H_
