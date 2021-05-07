// implementation file for SquareMatrix.h

// Vector

template<typename R, size_t n>
bool Vector<R,n>::operator==(const Vector<R,n> & other) const
{
  for (size_t i = 0; i < n; i++)
    if (this->v[i] != other[i]) return false;
  return true;
}

// for comparison in standard containers such as std::set and std::map
template<typename R, size_t n>
bool Vector<R,n>::operator<(const Vector<R,n> & other) const
{
  for (size_t i = 0; i < n; i++) {
    if (this->v[i] > other[i]) return false;
    if (this->v[i] < other[i]) return true;
  }
  return false;
}

template<typename R, size_t n>
Vector<R, n> Vector<R, n>::operator+(const Vector<R, n> & other) const {
  Vector<R, n> sum;
  for (size_t i = 0; i < n; i++)
    sum[i] = this->v[i] + other[i];
  return sum;
}

template<typename R, size_t n>
Vector<R, n> Vector<R, n>::operator*(const SquareMatrix<R, n>& mat) const
{
  Vector<R, n> prod;
  for (size_t i = 0; i < n; i++) {
    prod[i] = Math<R>::zero();
    for (size_t j = 0; j < n; j++)
      prod[i] += mat(j, i) * this->v[j];
  }
  return prod;
}

template<typename R, size_t n>
R Vector<R, n>::inner_product(const Vector<R, n> & vec1,
			      const Vector<R, n> & vec2)
{
  R prod = Math<R>::zero();

  for (size_t i = 0; i < n; i++)
    prod += vec1[i]*vec2[i];
  return prod;
}

template<typename R, typename S, size_t n>
FpElement<R, S> VectorFp<R, S, n>::inner_product(const VectorFp<R, S, n> & vec1,
						 const VectorFp<R, S, n> & vec2)
{
  FpElement<R, S> prod = Math<FpElement<R,S>>::zero();
  prod.set_field(vec1[0].field());

  for (size_t i = 0; i < n; i++)
    prod += vec1[i]*vec2[i];
  return prod;
}

// printing
template<typename R, size_t n>
std::ostream& operator<<(std::ostream& os, const Vector<R, n>& vec)
{
  os << "Vector(";
  for (size_t i = 0; i < n-1; i++)
    os << vec[i] << ",";
  os << vec[n-1] <<  ")";
  return os;
}

// SquareMatrix

template<typename R, size_t n>
void SquareMatrix<R, n>::deep_copy(const R mat[n][n])
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->mat[i][j] = mat[i][j];
}

// c-tors
template<typename R, size_t n>
SquareMatrix<R, n>::SquareMatrix(const R mat[n][n])
{ deep_copy(mat); }

template<typename R, size_t n>
SquareMatrix<R, n>::SquareMatrix(const SquareMatrix<R, n> & other)
{ deep_copy(other.mat); }

// assignment
template<typename R, size_t n>
SquareMatrix<R,n> &
SquareMatrix<R, n>::operator=(const SquareMatrix<R, n> & other)
{
  if (this !=  &other) {
    deep_copy(other.mat);
  }
  return (*this);
}

// arithmetic

template<typename R, size_t n>
SquareMatrix<R, n>
SquareMatrix<R, n>::operator*(const SquareMatrix<R, n>& other) const
{
  SquareMatrix<R, n> prod;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      prod(i,j) = 0;
      for (size_t k = 0; k < n; k++)
	prod(i,j) += this->mat[i][k]*other(k,j);
    }
  return prod;
}

template<typename R, size_t n>
Vector<R, n> SquareMatrix<R, n>::operator*(const Vector<R, n>& vec) const
{
  Vector<R, n> prod;
  for (size_t i = 0; i < n; i++) {
    prod[i] = Math<R>::zero();
    for (size_t j = 0; j < n; j++)
      prod[i] += this->mat[i][j] * vec[j];
  }
  return prod;
}

template<typename R, typename S, size_t n>
VectorFp<R, S, n>
SquareMatrixFp<R, S, n>::operator*(const VectorFp<R, S, n>& vec) const
{
  VectorFp<R, S, n> prod(this->GF);
  for (size_t i = 0; i < n; i++) {
    prod[i] = Math<R>::zero();
    for (size_t j = 0; j < n; j++)
      prod[i] += this->mat[i][j] * vec[j];
  }
  return prod;
}

template<typename R, size_t n>
SquareMatrix<R, n> SquareMatrix<R, n>::operator*(const R & scalar) const {
  SquareMatrix<R,n> prod;
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      prod(row,col) = scalar*this->mat[row][col];
  return prod;
}

template<typename R, size_t n>
SquareMatrix<R, n>  SquareMatrix<R, n>::operator/(const R & scalar) const {
  SquareMatrix<R, n> quo;
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      quo(row,col) = this->mat[row][col] / scalar;
  return quo;
}

// booleans
template<typename R, size_t n>
bool SquareMatrix<R, n>::operator==(const SquareMatrix<R, n>& other) const
{
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      if (mat[row][col] != other(row, col)) return false;
  return true;
}

// !! - TODO - replace == and < by compare to save time
template<typename R, size_t n>
bool SquareMatrix<R, n>::operator<(const SquareMatrix<R, n>& other) const
{
  for (size_t i = 0; i < n; i++) {
    if (mat[i][i] < other(i,i)) return true;
    if (mat[i][i] > other(i,i)) return false;
  }
  for (size_t col = 0; col < n-1; col++)
    for (size_t row = 0; row < n-1-col; row++) {
      if (abs(mat[row][row+col+1]) > abs(other(row, row+col+1))) return true;
      if (abs(mat[row][row+col+1]) < abs(other(row, row+col+1))) return false;
    }
  for (size_t col = 0; col < n-1; col++)
    for (size_t row = 0; row < n-1-col; row++) {
      if (mat[row][row+col+1] > other(row, row+col+1)) return true;
      if (mat[row][row+col+1] < other(row, row+col+1)) return false;
    }
  return false;
}

template<typename R, size_t n>
bool SquareMatrix<R, n>::is_upper_triangular() const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < i; j++)
      if (mat[i][j] != Math<R>::zero()) return false;
  return true;
}

template<typename R, size_t n>
bool SquareMatrix<R, n>::is_lower_triangular() const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = i+1; j < n; j++)
      if (mat[i][j] != Math<R>::zero()) return false;
  return true;
}

template<typename R, size_t n>
bool SquareMatrix<R, n>::is_symmetric() const
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < i; j++)
      if (mat[i][j] != mat[j][i]) return false;
  return true;
}

template<typename R, size_t n>
bool SquareMatrix<R, n>::is_positive_definite() const
{
  SquareMatrix<R, n> L;
  Vector<R, n> D;
  return cholesky(L,D);
}
  
// basic operations
template<typename R, size_t n>
SquareMatrix<R, n> SquareMatrix<R, n>::transpose(void) const
{
  SquareMatrix<R, n> trans;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      trans(i,j) = this->mat[j][i];
  return trans;
}

template<typename R, size_t n>
R SquareMatrix<R, n>::determinant(void) const
{
  // Instead of the previous ad-hoc method, we use Bareiss algorithm
  // to compute the determinant.
  // TODO - can do Cholesky, will be faster
  // !! TODO - Leibniz should be better when n <= 5 !?
  
  SquareMatrix<R, n+1> M;
  M(0,0) = 1;
  // init
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      M(row+1, col+1) = this->mat[row][col];
  for (size_t k = 1; k < n; k++)
    for (size_t i = k+1; i <= n; i++)
      for (size_t j = k+1; j <= n; j++) {
	if (M(k-1,k-1) == Math<R>::zero()) return Math<R>::zero();
	M(i,j) = (M(i,j)*M(k,k) - M(i,k)*M(k,j))/M(k-1,k-1);
      }
  return M(n,n);
}

// this solves using forward substitution in O(n^2)
// for lower triangular matrices
template<typename R, size_t n>
Vector<R, n>
SquareMatrix<R, n>::forward_substitution(const Vector<R,n> & vec) const
{
  Vector <R,n> sol;
  R sum;
  assert(is_lower_triangular());
 
  for (size_t i = 0; i < n; i++) {
    sum = Math<R>::zero();
    for (size_t j = 0; j < i; j++)
      sum += mat[i][j] * sol[j];
    sol[i] = (vec[i] - sum) / mat[i][i];
  } 
 
  return sol;
}

// this solves using forward substitution in O(n^2)
// for lower triangular matrices
template<typename R, size_t n>
SquareMatrix<R, n>
SquareMatrix<R, n>::inverse_lower_triangular(void) const
{
  SquareMatrix<R,n> inv = SquareMatrix<R, n>::identity();
  R sum;
  assert(is_lower_triangular());

  for (size_t col = 0; col < n; col++) {
    for (size_t i = col; i < n; i++) {
      sum = Math<R>::zero();
      for (size_t j = 0; j < i; j++)
	sum += mat[i][j] * inv(j, col);
      R delta = (i == col) ? Math<R>::one() : Math<R>::zero();
      inv(i,col) = (delta - sum) / mat[i][i];
    } 
  }
  return inv;
}

template<typename R, size_t n>
SquareMatrix<R, n>
SquareMatrix<R, n>::inverse_upper_triangular(void) const
{
  SquareMatrix<R,n> inv = SquareMatrix<R, n>::identity();
  R sum;
  assert(is_upper_triangular());

  for (size_t col = 0; col < n; col++) {
    for (size_t i = col+1; i > 0; i--) {
      sum = Math<R>::zero();
      for (size_t j = i; j < n; j++)
	sum += mat[i-1][j] * inv(j,col);
      R delta = (i-1 == col) ? Math<R>::one() : Math<R>::zero();
      inv(i-1,col) = (delta - sum) / mat[i-1][i-1];
    } 
  }
  return inv;
}

// this solves using backward substitution in O(n^2)
// for upper triangular matrices
template<typename R, size_t n>
Vector<R, n>
SquareMatrix<R, n>::backward_substitution(const Vector<R,n> & vec) const
{
  Vector <R,n> sol;
  R sum;
  assert(is_upper_triangular());
  
  for (size_t i = n; i > 0; i--) {
    sum = Math<R>::zero();
    for (size_t j = i; j < n; j++)
      sum += mat[i-1][j] * sol[j];
    sol[i-1] = (vec[i-1] - sum) / mat[i-1][i-1];
  } 
  
  return sol;
}

// returns false if the matrix is not positive definite
template<typename R, size_t n>
bool
SquareMatrix<R, n>::cholesky(SquareMatrix<R, n>& L,  Vector<R,n> & D) const
{
  assert(is_symmetric());
  L = SquareMatrix<R, n>::identity();
  R sum;
  for (size_t j = 0; j < n; j++) {
    sum = Math<R>::zero();
    for (size_t k = 0; k < j; k++)
      sum += L(j,k)*L(j,k)*D[k];
    D[j] = mat[j][j] - sum;
    if (D[j] == Math<R>::zero()) return false;
    for (size_t i = j+1; i < n; i++) {
      sum = Math<R>::zero();
      for (size_t k = 0; k < j; k++)
	sum += L(i,k)*L(j,k)*D[k];
      L(i,j) = (mat[i][j] - sum) / D[j];
    } 
  }
  return true;
}

// This solves only for symmetric positive definite matrices
// using LDL
template<typename R, size_t n>
Vector<R, n>
SquareMatrix<R, n>::solve(const Vector<R,n> & vec) const
{
  assert(is_symmetric());
  Vector<R, n> sol;
  SquareMatrix<R, n> L;
  Vector<R, n> D;
  bool is_positive_definite = cholesky(L, D);
  assert(is_positive_definite);
  sol = L.forward_substitution(vec);
  for (size_t i = 0; i < n; i++)
    sol[i] /= D[i];
  sol = L.transpose().backward_substitution(sol);
  return sol;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::swap_rows(size_t row1, size_t row2)
{
  Vector<R, n> temp_row;

  for (size_t col = 0; col < n; col++)
    temp_row[col] = mat[row1][col];

  for (size_t col = 0; col < n; col++)
    mat[row1][col] = mat[row2][col];

  for (size_t col = 0; col < n; col++)
    mat[row2][col] = temp_row[col];

  return;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::swap_cols(size_t col1, size_t col2)
{
  Vector<R, n> temp_col;

  for (size_t row = 0; row < n; row++)
    temp_col[row] = mat[row][col1];

  for (size_t row = 0; row < n; row++)
    mat[row][col1] = mat[row][col2];

  for (size_t row = 0; row < n; row++)
    mat[row][col2] = temp_col[row];

  return;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::multiply_row(size_t row, const R & val)
{
  for (size_t col = 0; col < n; col++)
    mat[row][col] *= val;
  return;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::multiply_col(size_t col, const R & val)
{
  for (size_t row = 0; row < n; row++)
    mat[row][col] *= val;
  return;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::add_row(size_t row_to, size_t row_from, const R & val)
{
  for (size_t col = 0; col < n; col++) {
    mat[row_to][col] += val * mat[row_from][col];
  }
  return;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::add_col(size_t col_to, size_t col_from, const R & val)
{
  for (size_t row = 0; row < n; row++) {
    mat[row][col_to] += val * mat[row][col_from];
  }
  return;
}
  
// a general one, just in case
template<typename R, size_t n>
SquareMatrix<R, n> SquareMatrix<R, n>::inverse(void) const
{
  if (is_lower_triangular()) return inverse_lower_triangular();
  if (is_upper_triangular()) return inverse_upper_triangular();
  if (is_symmetric()) {
    SquareMatrix<R, n> L;
    Vector<R, n> D;
    bool is_positive_definite = cholesky(L,D);
    if (is_positive_definite) {
      SquareMatrix<R,n> L_inv = L.inverse();
      SquareMatrix<R,n> L_inv_t = L_inv.transpose();
      for (size_t i = 0; i < n; i++)
	for (size_t j = 0; j < n; j++)
	  L_inv(i,j) /= D[i];
      
      return L_inv_t * L_inv;
    }
  }
  SquareMatrix<R, n> inv = identity();
  SquareMatrix<R, n> echelon(mat);
  size_t pivot_row = 0;
  size_t pivot_col = 0;
  size_t row_max;
  R max_val, val, factor;
  
  while ((pivot_row < n) && (pivot_col < n)) {
    row_max = pivot_row;
    max_val = abs(echelon(row_max, pivot_col));
    for (size_t row = pivot_row+1; row < n; row++) {
      val = abs(echelon(row, pivot_col));
      if (max_val < val) {
	row_max = row;
	max_val = val;
      }
    }
    if (max_val == 0) {
      pivot_col++;
    }
    else {
      echelon.swap_rows(pivot_row, row_max);
      inv.swap_rows(pivot_row, row_max);
      R scalar = 1 / echelon(pivot_row, pivot_col);
      echelon.multiply_row(pivot_row, scalar);
      inv.multiply_row(pivot_row, scalar);
      // for reduced row echelon form we need also the rows before
      for (size_t row = 0; row < pivot_row; row++) {
	factor = echelon(row, pivot_col);
	echelon(row, pivot_col) = 0;
	for (size_t col = pivot_col + 1; col < n; col++) {
	  echelon(row,col) -= factor * echelon(pivot_row, col);
	}
	for (size_t col = 0; col < n; col++) {
	  inv(row, col) -= factor * inv(pivot_row, col);
	}
      }
      for (size_t row = pivot_row+1; row < n; row++) {
	// factor = echelon(row,pivot_col) / echelon(pivot_row, pivot_col);
	factor = echelon(row, pivot_col);
	echelon(row, pivot_col) = 0;
	for (size_t col = pivot_col + 1; col < n; col++) {
	  echelon(row,col) -= factor * echelon(pivot_row, col);
	}
	for (size_t col = 0; col < n; col++) {
	  inv(row, col) -= factor * inv(pivot_row, col);
	}
      }
      
      pivot_row++;
      pivot_col++;
    }
  }
#ifdef DEBUG
  assert((*this)*inv == identity());
#endif
  return inv;
}

// static functions
template<typename R, size_t n>
Rational<R>
SquareMatrix<R,n>::inner_product(const SquareMatrix<R, n> & F,
				 const SquareMatrix<Rational<R>, n> & S,
				 size_t idx1, size_t idx2)
{
  Rational<R> ans = Math< Rational<R> >::zero();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ans += S(idx1, i) * F(i,j) * S(idx2, j);
  return ans;
}

// global constants
template<typename R, size_t n>
SquareMatrix<R, n> SquareMatrix<R, n>::identity(void)
{
  SquareMatrix<R, n> id;
  R one = Math<R>::one();
  R zero = Math<R>::zero();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      id(i,j) = (i == j) ? one : zero;
  return id; 
}

// printing
template<typename R, size_t n>
std::ostream& operator<<(std::ostream& os, const SquareMatrix<R, n>& a)
{
  os << "Matrix(Integers(), " << n << " , (";
  for (size_t row = 0; row < n-1; row++) {
    for (size_t col = 0; col < n; col++)
      os << a(row, col) << ",";
  }
  for (size_t col = 0; col < n-1; col++)
    os << a(n-1,col) << ",";
  os << a(n-1,n-1) <<  "))";
    
  return os;
}

template<typename R, size_t n>
std::ostream & Vector<R,n>::pretty_print(std::ostream & os,
					 size_t upTo) const
{
  for (size_t i = 0; i < upTo; i++) {
    os << (*this)[i] << " ";
  }
  os << std::endl;
 
  return os;
}

template<typename R, size_t n>
std::ostream & SquareMatrix<R,n>::pretty_print(std::ostream & os,
					       size_t upTo) const
{
 for (size_t row = 0; row < upTo-1; row++) {
    for (size_t col = 0; col < upTo; col++)
      os << (*this)(row, col) << " ";
    os << std::endl;
  }
  for (size_t col = 0; col < upTo-1; col++)
    os << (*this)(upTo-1,col) << " ";
  os << (*this)(upTo-1,upTo-1) <<  std::endl;
    
  return os;
}
