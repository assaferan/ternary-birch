#ifndef __MATRIX_H_
#define __MATRIX_H_

#include <vector>

#include "Polynomial.h"

template<typename R>
class Matrix
{
public:
  Matrix(const std::vector<R> & data, size_t nrows, size_t ncols)
    : nrows_(nrows), ncols_(ncols), data_(data) {}
  template <size_t n>
  Matrix(const R data[n][n]);
  template <size_t n>
  Matrix(const SquareMatrix<R, n> &);
  Matrix(size_t nrows, size_t ncols)
    : nrows_(nrows), ncols_(ncols), data_(nrows*ncols) {}
  const R & operator()(size_t row, size_t col) const
  {
#ifdef DEBUG
    assert( (ncols_*row+col < nrows_*ncols_) && (0 <= ncols_*row+col) );
#endif    
    return data_[ncols_*row+col];
  }
  R & operator()(size_t row, size_t col)
  {
#ifdef DEBUG
    assert( (ncols_*row+col < nrows_*ncols_) && (0 <= ncols_*row+col) );
#endif
    return data_[ncols_*row+col];
  }
  size_t nrows() const {return nrows_;}
  size_t ncols() const {return ncols_;}

  // return the i-th row
  std::vector<R> operator[](size_t i) const;
  
  R determinant() const;

  size_t rank() const;

  Matrix<R> kernel() const;

  Matrix<R> left_kernel() const;

  // restrict matrix to the subspace specified by the argument
  Matrix<R> restrict(const Matrix<R> & ) const;

  R trace() const;
  
  UnivariatePoly<Z> char_poly() const;
  
  static Matrix<R> diagonal_join(const std::vector< Matrix<R> > & mats);

  static Matrix<R> identity(size_t n);
  
  // TODO - just change access resolution to the same vector instead
  Matrix<R> transpose() const;

  void swap_rows(size_t, size_t);

  Matrix<R> operator+(const Matrix<R> &) const;
  Matrix<R> operator-(const Matrix<R> &) const;
  Matrix<R> operator*(const Matrix<R> &) const;

  Matrix<R>& operator+=(const Matrix<R> &);
  Matrix<R>& operator-=(const Matrix<R> &);
  Matrix<R>& operator*=(const Matrix<R> &);

  Matrix<R> operator*(const R & a) const;

  // in-place row-echelon form for the matrix echelon,
  // returns the rank and the transformation matrix trans
  static size_t row_echelon(Matrix<R> & echelon, Matrix<R>& trans);
  
protected:
  size_t nrows_;
  size_t ncols_;
  std::vector<R> data_;
};

template<typename R>
Matrix<R> operator*(const R & a, const Matrix<R> & mat)
{ return mat*a; }

template<typename R, typename S>
class MatrixFp : public Matrix<FpElement<R, S> >
{
public:

  MatrixFp(std::shared_ptr<const Fp<R,S>> GF, size_t nrows, size_t ncols)
    : Matrix<FpElement<R,S> >(nrows, ncols)
  { set_field(GF);}

  MatrixFp(std::shared_ptr<const Fp<R,S>> GF,
	   const Matrix<FpElement<R,S> > & other)
    : Matrix<FpElement<R,S> >(other)
  { set_field(GF);}

  size_t rank() const;
  MatrixFp<R, S> kernel() const;
  MatrixFp<R, S> left_kernel() const;
  
protected:
  std::shared_ptr<const Fp<R,S>> GF;

  void set_field(std::shared_ptr<const Fp<R,S>> GF)
  {
    this->GF = GF;
    size_t idx = 0;
    for (size_t i = 0; i < this->nrows(); i++)
      for (size_t j = 0; j < this->ncols(); j++)
	this->data_[idx++].set_field(GF);
  }
};

template <typename R>
std::ostream& operator<<(std::ostream & os, const Matrix<R> & mat)
{
  for (size_t i = 0; i < mat.nrows(); i++) {
    for (size_t j = 0; j < mat.ncols(); j++)
      os << mat(i,j) << " ";
    os << std::endl;
  }
  return os;
}

#include "Matrix.inl"

#endif // __MATRIX_H_
