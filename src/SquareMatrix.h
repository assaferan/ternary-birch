#ifndef __SQUAREMATRIX_H_
#define __SQUAREMATRIX_H_

#include "Rational.h"

template<typename R, size_t n>
class Vector {
public:

  // c-tor - default constructor constructs the zero vector
  Vector()
  {
    for (size_t i = 0; i < n; i++) this->v[i] = Math<R>::zero();
  }
  // access
  const R& operator[](size_t i) const {return v[i]; }
  R& operator[](size_t i) {return v[i];}

  // arithmetic
  Vector<R, n> operator-() const;
  Vector<R, n> operator+(const Vector<R, n> &) const;
  Vector<R, n> operator-(const Vector<R, n> &) const;
  Vector<R, n> operator*(const R & a) const;

  Vector<R, n> & operator+=(const Vector<R,n> &);
  Vector<R, n> & operator-=(const Vector<R,n> &);
  
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

  bool operator!=(const Vector<R,n> & other) const
  {return !((*this)==other);}

  bool operator<(const Vector<R,n> &) const;

  bool is_zero(void) const
  {
    for (size_t i = 0; i < n; i++)
      if (v[i] != 0)
	return false;
    return true;
  }

  std::ostream & pretty_print(std::ostream &, size_t upTo = n) const;
  
protected:
  R v[n];
};

//scalar multiple
template<typename R, size_t n>
Vector<R,n> operator*(const R & a, const Vector<R,n> & v)
{ return v*a; }

template<typename R, typename S, size_t n>
class VectorFp : public Vector<FpElement<R, S>, n>
{
public:
  VectorFp(std::shared_ptr<const Fp<R,S>> GF)
  {
    this->GF = GF;
    for (size_t i = 0; i < n; i++)
      this->v[i].set_field(GF);
  }

  static FpElement<R, S> inner_product(const VectorFp<R, S, n> &,
				       const VectorFp<R, S, n> &);
protected:
  std::shared_ptr<const Fp<R,S>> GF;
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

  // return the i-th row
  Vector<R, n> operator[](size_t i) const;
  
  // arithmetic
  SquareMatrix<R, n> operator*(const SquareMatrix<R, n>&) const;
  Vector<R, n> operator*(const Vector<R, n>& vec) const;
  SquareMatrix<R, n> operator*(const R &) const;
  SquareMatrix<R, n> operator/(const R &) const;

  // booleans
  bool operator==(const SquareMatrix<R, n>&) const;
  bool operator!=(const SquareMatrix<R, n>& other) const
  { return !((*this) == other);}
  // ordering of the matrices for Minkowski reduction
  bool operator<(const SquareMatrix<R, n>&) const;
  bool is_upper_triangular() const;
  bool is_lower_triangular() const;
  bool is_symmetric() const;
  bool is_positive_definite() const;
  
  // basic operations
  void set_identity(void);
  
  SquareMatrix<R, n> transpose(void) const;
  // !! TODO - save the inverse and track it
  // to save computation
  SquareMatrix<R, n> inverse(void) const;
  template<size_t m>
  SquareMatrix<R, m> submatrix(size_t idxs[m]) const;
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
  static Rational<R> inner_product(const SquareMatrix<R,n> & F,
				   const SquareMatrix<Rational<R>,n> & S,
				   size_t idx1, size_t idx2);
  
  // global constants
  static SquareMatrix<R, n> identity(void);

  std::ostream & pretty_print(std::ostream &, size_t upTo = n) const;
  
protected:
  R mat[n][n];
  
  // helper functions
  void deep_copy(const R mat[n][n]);
  
  Vector<R, n> forward_substitution(const Vector<R,n> & vec) const;
  Vector<R, n> backward_substitution(const Vector<R,n> & vec) const;
  SquareMatrix<R, n> inverse_lower_triangular(void) const;
  SquareMatrix<R, n> inverse_upper_triangular(void) const;
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
	      fnv = (fnv ^ vec[i]) * FNV_PRIME;
            
            return fnv;
        }
    };
}

template<typename R, typename S, size_t n>
class SquareMatrixFp : public SquareMatrix<FpElement<R, S>, n>
{
public:
  SquareMatrixFp(std::shared_ptr<const Fp<R,S>> GF)
  { set_field(GF);}

  SquareMatrixFp(std::shared_ptr<const Fp<R,S>> GF,
		 const SquareMatrix< FpElement<R,S>, n> & other)
    : SquareMatrix< FpElement<R,S>, n>(other)
  { set_field(GF); }

  // return the i-th row
  VectorFp<R, S, n> operator[](size_t i) const;
  
  VectorFp<R, S, n> operator*(const VectorFp<R, S, n>& ) const;
  SquareMatrixFp<R, S, n> operator*(const SquareMatrixFp<R, S, n>&) const;
  
  std::shared_ptr<const Fp<R,S>> field(void) const
  {return GF;}

  SquareMatrixFp<R, S, n> transpose(void) const
  {
    SquareMatrix< FpElement<R, S>, n> tr =
      SquareMatrix< FpElement<R, S>, n>::transpose();
    return SquareMatrixFp<R, S, n>(this->GF, tr);
  }

  void set_identity(void);
  
protected:
  std::shared_ptr<const Fp<R,S>> GF;

  void set_field(std::shared_ptr<const Fp<R,S>> GF)
  {
    this->GF = GF;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	this->mat[i][j].set_field(GF);
  }
};

#include "SquareMatrix.inl"

#endif // __SQUAREMATRIX_H_
