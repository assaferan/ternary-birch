#ifndef __MATH_H_
#define __MATH_H_

#include "birch.h"

template<typename R>
class Rational;

template<typename R>
class Math
{
public:
  static int hilbert_symbol(R a, R b, const R& p);
  static int kronecker_symbol(const R & a, const R & n);
  static size_t valuation(const R& a, const R& p);
  static bool is_local_square(const R& a, const R& p);
  static std::vector< std::pair<R, size_t> > factorization(const R & num);
  static bool is_square(const R & num);
  static Rational<R> bernoulli_number(const size_t &);
  static R binomial_coefficient(const R & n, const R & k);
  static std::vector< Rational<R> > bernoulli_poly(const size_t & n);
  static Rational<R> bernoulli_number(const size_t & n, const R & d);
};

template<typename R>
class Matrix
{
public:
  Matrix(const std::vector<R> & data, size_t nrows, size_t ncols)
    : nrows_(nrows), ncols_(ncols), data_(data) {}
  Matrix(const R** data, size_t nrows, size_t ncols)
    : nrows_(nrows), ncols_(ncols), data_(nrows*ncols)
  { idx := 0;
    for (size_t row = 0; row < nrows; row++)
      for (size_t col = 0; col < ncols; col++)
	data_[idx++] = data[row][col];
  }
  Matrix(size_t nrows, size_t ncols)
    : nrows_(nrows), ncols_(ncols), data_(nrows*ncols) {}
  const R & operator()(size_t row, size_t col) const
  {return data_[_ncols*row+col];}
  R & operator()(size_t row, size_t col)
  {return data_[_ncols*row+col];}
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
        M[row+1][col+1] = (*this)[row][col];
    for (size_t k = 1; k < n; k++)
      for (size_t i = k+1; i <= n; i++)
        for (size_t j = k+1; j <= n; j++)
          M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j])/M[k-1][k-1];
    return M[n][n];
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
	  diag[big_row + row][big_col + col] = mat[row][col];
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
	data_t[idx++] = (*this)[col][row];
    Matrix<R> tr(data_t, ncols_, nrows_);
    return tr;
  }
  Matrix<R> operator*(const Matrix<R> & other) const
  {
    size_t nrows = this->nrows_;
    size_t ncols = other.ncols_;
    std::vector<R> data(nrows*ncols);
    size_t idx = 0;
    assert( this->ncols_ == other.nrows_ );
    for (size_t row = 0; row < nrows_; row++)
      for (size_t col = 0; col < ncols_; col++) {
	data_t[idx] = 0;
	for (size_t j = 0; j < this->ncols_; j++)
	  data_t[idx] += (*this)[row][j]*(*this)[j][col];
	idx++;
      }
    Matrix<R> prod(data_t, nrows, ncols);
    return prod;
  }
protected:
  size_t nrows_;
  size_t ncols_;
  std::vector<R> data_;
};

template<typename R>
class Rational
{
public:
  // c-tors
  Rational(const R& num, const R& denom) : num_(num), denom_(denom) {}

  Rational(const R[2] & nums) : num_(nums[0]), denom_(nums[1]) {}

  Rational(const R & num) : num_(num), denom_(1) {}
  
  // default c-tor
  Rational() : num_(0), denom_(1) {}

  // arithmetic
  Rational<R> operator+() const {return Rational(num_, denom_); }
  Rational<R> operator-() const {return Rational(-num_, denom); }
  Rational<R> operator+(const Rational<R> &) const;
  Rational<R> operator-(const Rational<R> &b) const {return (*this)+(-b); }
  Rational<R> operator*(const Rational<R> &) const;
  Rational<R> operator/(const Rational<R> &) const;

  // comparison
  bool operator==(const Rational<R> &) const;
  bool operator!=(const Rational<R> &b) const {return !((*this)==b); }
  bool operator<(const Rational<R> &) const;
  bool operator>(const Rational<R> &b) const {return b < (*this); }
  bool operator<=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  bool operator>=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) > b); }
  
protected:
  R num_;
  R denom_;

private:
  void reduce(void);
};

#endif // __MATH_H_
