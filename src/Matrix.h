#ifndef __MATRIX_H_
#define __MATRIX_H_

#include <vector>

template<typename R>
class Matrix
{
public:
  Matrix(const std::vector<R> & data, size_t nrows, size_t ncols)
    : nrows_(nrows), ncols_(ncols), data_(data) {}
  template <size_t n>
  Matrix(const R data[n][n])
    : nrows_(n), ncols_(n), data_(n*n)
  { size_t idx = 0;
    for (size_t row = 0; row < nrows_; row++)
      for (size_t col = 0; col < ncols_; col++)
	data_[idx++] = data[row][col];
  }
  Matrix(size_t nrows, size_t ncols)
    : nrows_(nrows), ncols_(ncols), data_(nrows*ncols) {}
  const R & operator()(size_t row, size_t col) const
  {return data_[ncols_*row+col];}
  R & operator()(size_t row, size_t col)
  {return data_[ncols_*row+col];}
  size_t nrows() const {return nrows_;}
  size_t ncols() const {return ncols_;}

  R determinant() const
  {
    assert(nrows_ == ncols_);
    size_t n = nrows_;
    Matrix<R> M(n+1, n+1);
    M(0,0) = 1;
    // init
    for (size_t row = 0; row < n; row++)
      for (size_t col = 0; col < n; col++)
        M(row+1, col+1) = (*this)(row, col);
    for (size_t k = 1; k < n; k++)
      for (size_t i = k+1; i <= n; i++)
        for (size_t j = k+1; j <= n; j++)
          M(i,j) = (M(i,j)*M(k,k) - M(i,k)*M(k,j))/M(k-1,k-1);
    return M(n,n);
  }
  
  static Matrix<R> diagonal_join(const std::vector< Matrix<R> > & mats)
  {
    size_t nrows = 0;
    size_t ncols = 0;
    for (Matrix<R> mat : mats) {
      nrows += mat.nrows();
      ncols += mat.ncols();
    }
    Matrix<R> diag(nrows, ncols);
    size_t big_row = 0;
    size_t big_col = 0;
    for (Matrix<R> mat : mats) {
      for (size_t row = 0; row < mat.nrows(); row++)
	for (size_t col = 0; col < mat.ncols(); col++)
	  diag(big_row + row, big_col + col) = mat(row, col);
      big_row += mat.nrows();
      big_col += mat.ncols();
    }
    return diag;
  }
  
  // TODO - just change access resolution to the same vector instead
  Matrix<R> transpose() const
  {
    std::vector<R> data_t(nrows_*ncols_);
    size_t idx = 0;
    for (size_t row = 0; row < nrows_; row++)
      for (size_t col = 0; col < ncols_; col++)
	data_t[idx++] = (*this)(col, row);
    Matrix<R> tr(data_t, ncols_, nrows_);
    return tr;
  }
  Matrix<R> operator*(const Matrix<R> & other) const
  {
    size_t nrows = this->nrows_;
    size_t ncols = other.ncols_;
    assert( this->ncols_ == other.nrows_ );
    Matrix<R> prod(nrows, ncols);
    
    for (size_t row = 0; row < nrows_; row++)
      for (size_t col = 0; col < ncols_; col++) {
	prod(row,col) = 0;
	for (size_t j = 0; j < this->ncols_; j++)
	  prod(row,col) += (*this)(row,j)*(*this)(j,col);
      }
    return prod;
  }
protected:
  size_t nrows_;
  size_t ncols_;
  std::vector<R> data_;
};

template <typename R>
std::ostream& operator<<(std::ostream os&, const Matrix<R> & mat)
{
  for (size_t i = 0; i < mat.nrows(); i++) {
    for (size_t j = 0; j < mat.ncols(); j++)
      os << mat[i][j] << " ";
    os << std::endl;
  }
}

#include "Matrix.inl"

#endif // __MATRIX_H_
