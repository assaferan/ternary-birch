#ifndef __SQUAREMATRIX_H_
#define __SQUAREMATRIX_H_

#include "Rational.h"

template<typename R, size_t n>
class Vector {
public:

  // access
  const R& operator[](size_t i) const {return v[i]; }
  R& operator[](size_t i) {return v[i];}

  // arithmetic
  Vector<R, n> operator+(const Vector<R, n> &) const;
  
  // considering the vector as a row vector
  Vector<R, n> operator*(const SquareMatrix<R, n>& mat) const;

  // inner product
  static R inner_product(const Vector<R, n> & , const Vector<R, n> &);

  // assignment
  Vector<R, n> & operator=(const Vector<R, n>& other)
  { if (this != (&other)) {
      for (size_t i = 0; i < n; i++) this->v[i] = other[i];
    }
    return (*this);
  }
  
  // booleans
  bool operator==(const Vector<R,n> &) const;
  
protected:
  R v[n];
};

// printing
template<typename R, size_t n>
std::ostream& operator<<(std::ostream&, const Vector<R, n>&);

template<typename R, size_t n>
class SquareMatrix {
public:
  // c-tor
  SquareMatrix() = default;

  SquareMatrix(const R mat[n][n]);
  SquareMatrix(const SquareMatrix<R, n> & mat);

  // assignment
  SquareMatrix<R,n> & operator=(const SquareMatrix<R,n> &);
  
  // access
  const R& operator()(size_t i, size_t j) const {return mat[i][j]; }
  R& operator()(size_t i, size_t j) {return mat[i][j];}

  // arithmetic
  SquareMatrix<R, n> operator*(const SquareMatrix<R, n>&) const;
  Vector<R, n> operator*(const Vector<R, n>& vec) const;
  SquareMatrix<R, n> operator*(const R &) const;
  SquareMatrix<R, n> operator/(const R &) const;

  // booleans
  bool operator==(const SquareMatrix<R, n>&) const;
  // ordering of the matrices for Minkowski reduction
  bool operator<(const SquareMatrix<R, n>&) const;
  bool is_upper_triangular() const;
  bool is_lower_triangular() const;
  bool is_symmetric() const;
  bool is_positive_definite() const;
  
  // basic operations
  SquareMatrix<R, n> transpose(void) const;
  // !! TODO - save the inverse and track it
  // to save computation
  SquareMatrix<R, n> inverse(void) const;
  R determinant(void) const;

  Vector<R,n> solve(const Vector<R,n> & vec) const;

  // more complex operations that might be useful outside the class
  bool cholesky(SquareMatrix<R, n>& L,  Vector<R,n> & D) const;
  
  // elementary operations
  void swap_rows(size_t row1, size_t row2);
  void swap_cols(size_t col1, size_t col2);
  void multiply_row(size_t row, const R & val);
  void multiply_col(size_t col, const R & val);
  void add_row(size_t row_to, size_t row_from, const R & val);
  void add_col(size_t col_to, size_t col_from, const R & val);
    
  // static functions
  // compute S[idx1]*F*S[idx2]^t
  // this is needed in this form for jordan decomposition
  Rational<R> inner_product(const SquareMatrix<R,n> & F,
			    const SquareMatrix<Rational<R>,n> & S,
			    size_t idx1, size_t idx2);
  
  // global constants
  static SquareMatrix<R, n> identity(void);
  
protected:
  R mat[n][n];
  
  // helper functions
  void deep_copy(const R & mat[n][n]);
  
  Vector<R, n> forward_substitution(const Vector<R,n> & vec) const;
  Vector<R, n> backward_substitution(const Vector<R,n> & vec) const;
  
};

// printing
template<typename R, size_t n>
std::ostream& operator<<(std::ostream&, const SquareMatrix<R, n>&);

namespace std
{
  template<typename R, size_t n>
  struct hash<Vector<R, n> >
    {
      Z64 operator()(const Vector<R,n>& vec) const
        {
            Z64 fnv = FNV_OFFSET;
	    for (size_t i = 0; i < n; i++)
	      fnv = (fnv ^ vec.v[i]) * FNV_PRIME;
            
            return fnv;
        }
    };
}

#include "SquareMatrix.inl"

#endif // __SQUAREMATRIX_H_
