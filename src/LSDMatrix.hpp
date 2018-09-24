/// LSDMatrix
///
/// LSDMatrix is a basic, simple matrix class designed to store elements
/// contiguously in memory. No fancy tools (matrix algebra etc).
///
/// Experimental use!
///
/// @author DAV, 2016

#ifndef LSDMATRIX_HPP
#define LSDMATRIX_HPP


// http://www.graphics.cornell.edu/~martin/docs/c++-faq/freestore-mgmt.html#[16.15]
// Class Template version from the C++FAQ
// I think this is my favourite
/// @brief Simple, barebones, matrix class with dimensions that can be set at run time,
/// array is allocated contiguously in memory.
template<class T>
class LSDMatrix2D
{
public:
  LSDMatrix2D(unsigned rows, unsigned ncols);

  // For size being zero, we should throw error
  class BadSize { };

  // Law of big 3
  ~LSDMatrix2D();
  LSDMatrix2D(const LSDMatrix2D<T>& m);
  LSDMatrix2D& operator= (const LSDMatrix2D<T>& m);

  // Array access methods to get element with (i,j) notation
  T& operator() (unsigned i, unsigned j);
  const T& operator() (unsigned i, unsigned j) const;

  // Throw bounds violation if i, j too big
  class BoundsViolation { };

private:
  T* data_;
  unsigned nrows_, ncols_;
};

// Access matrix elements with m(i,j) notation
template<class T>
inline T& LSDMatrix2D<T>::operator() (unsigned row, unsigned col)
{
  if (row >= nrows_ || col >= ncols_) throw BoundsViolation();
  return data_[row*ncols_ + col];
}

// Access matrix elements with m(i,j) notation (constatnt)
template<class T>
inline const T& LSDMatrix2D<T>::operator() (unsigned row, unsigned col) const
{
  if (row >= nrows_ || col >= ncols_) throw BoundsViolation();
  return data_[row*ncols_ + col];
}

// Declare matrix
template<class T>
inline LSDMatrix2D<T>::LSDMatrix2D(unsigned rows, unsigned ncols)
  : data_ (new T[nrows * ncols]),
    nrows_ (nrows),
    ncols_ (ncols)
{
  if (nrows == 0 || ncols == 0)
    throw BadSize();
}

// Clean up after we're done with our matrix!
template<class T>
inline LSDMatrix2D<T>::~LSDMatrix2D()
{
  delete[] data_;
}


// OTHER EXAMPLES
// http://stackoverflow.com/a/28841507/1953517
class LSDArray2D
{
  int* array;
  int m_width;
public:
  LSDArray2D( int w, int h )
    :m_width( w ), array( new int[ w * h ] ) {}

  ~LSDArray2D() { delete[] array; }

  int at( int x, int y )
  const { return array[ index( x, y ) ]; }

protected:
  int index( int x, int y )
  const { return x + m_width * y; }

};


// http://stackoverflow.com/a/32279494/1953517
#include <memory>

class LSDGrid2D
{
    size_t _rows;
    size_t _columns;
    std::unique_ptr<int[]> data;

public:

    LSDGrid2D(size_t rows, size_t columns)
        : _rows{rows}, _columns{columns}, data{std::make_unique<int[]>(rows * columns)} {
    }


    size_t rows() const {
        return _rows;
    }

    size_t columns() const {
        return _columns;
    }

    int * operator[](size_t row) {
        return row * _columns + data.get();
    }

};




// Example from the C++ cookbook
// Require recipe 11.12!
#include <valarray>
#include <numeric>
#include <algorithm>

template<class Value_T>
class LSDMatrix
{
public:
  typedef Value_T value_type;
  typedef LSDMatrix self;
  typedef value_type* iterator;
  typedef const value_type* const_iterator;
  typedef Value_T* row_type;
  typedef stride_iter<value_type*> col_type; // need to implement stride iter!
  typedef const value_type* const_row_type;
  typedef stride_iter<const value_type*> const_col_type;

  // CONSTRUCTORS
  LSDMatrix() : nrows(0), ncols(0, m() { }

  LSDMatrix(int r, int c)
  : nrows(r), ncols(c), m(r*c) { }

  LSDMatrix(const self& x)
  : m(x.m), nrows(x.rows), ncols(x.cols) { }

  template<typename T>
  explicit LSDMatrix(const valarray<T> & x)
  : m(x.size() + 1), nrows(x.size()), ncols(1)
  {
    for(int i=0; i<x.size(); i++)
    {
      m[i] = x[i];
    }
  }

  // Allow construction from matrices of other types
  template<typename T>
  explicit LSDMatrix(const LSDMatrix<T>& x)
  : m(x.size() + 1), nrows(x.rows), ncols(x.cols)
  {
    copy(x.begin(), x.end(), m.begin());
  }

  // Public Functions
  int rows() const { return nrows; }
  int cols() const { return ncols; }
  int size() const { return nrows * ncols; }

  // Element access
  row_type row_begin(int n)
  {
    return &m[n * ncols];
  }

  row_type row_end


};

#endif // LSDMATRIX_HPP

