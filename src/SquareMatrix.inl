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
Vector<R, n> Vector<R, n>::operator-() const {
  Vector<R, n> neg;
  for (size_t i = 0; i < n; i++)
    neg[i] = -this->v[i];
  return neg;
}

template<typename R, size_t n>
Vector<R, n> Vector<R, n>::operator+(const Vector<R, n> & other) const {
  Vector<R, n> sum;
  for (size_t i = 0; i < n; i++)
    sum[i] = this->v[i] + other[i];
  return sum;
}

template<typename R, size_t n>
Vector<R, n> Vector<R, n>::operator-(const Vector<R, n> & other) const {
  Vector<R, n> diff;
  for (size_t i = 0; i < n; i++)
    diff[i] = this->v[i] - other[i];
  return diff;
}

template<typename R, size_t n>
Vector<R, n> Vector<R, n>::operator*(const R & a) const {
  Vector<R, n> prod;
  for (size_t i = 0; i < n; i++)
    prod[i] = a * this->v[i];
  return prod;
}

template<typename R, size_t n>
Vector<R, n> & Vector<R, n>::operator+=(const Vector<R,n> & other)
{
  for (size_t i = 0; i < n; i++)
    this->v[i] += other[i];
  return (*this);
}

template<typename R, size_t n>
Vector<R, n> & Vector<R, n>::operator-=(const Vector<R,n> & other)
{
  for (size_t i = 0; i < n; i++)
    this->v[i] -= other[i];
  return (*this);
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

// access
template<typename R, size_t n>
Vector<R, n> SquareMatrix<R, n>::operator[](size_t i) const
{
#ifdef DEBUG
  assert(i < n);
#endif 
  Vector<R, n> v;
  for (size_t j = 0; j < n; j++)
    v[j] = (*this)(i,j);
  return v;
}

template<typename R, typename S, size_t n>
VectorFp<R, S, n> SquareMatrixFp<R, S, n>::operator[](size_t i) const
{
#ifdef DEBUG
  assert(i < n);
#endif 
  VectorFp<R, S, n> v(this->GF);
  for (size_t j = 0; j < n; j++)
    v[j] = (*this)(i,j);
  return v;
}

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
      prod(i,j) = Math<R>::zero();
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
SquareMatrix<R, n>  SquareMatrix<R, n>::operator/(const R & scalar) const
{
#ifdef DEBUG
  assert(scalar != 0);
#endif 
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
  if (!is_symmetric()) {
    for (size_t col = 0; col < n-1; col++)
      for (size_t row = 0; row < n-1-col; row++) {
	if (mat[row+col+1][col] > other(row+col+1, col)) return true;
	if (mat[row+col+1][col] < other(row+col+1, col)) return false;
      }
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
  return ldl(L,D);
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
template<size_t m>
SquareMatrix<R, m> SquareMatrix<R, n>::submatrix(size_t idxs[m]) const
{
  SquareMatrix<R, m> sub;
  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < m; j++)
      sub(i,j) = (*this)(idxs[i], idxs[j]);
  
  return sub;
}

template<typename R, size_t n>
R SquareMatrix<R, n>::determinant(void) const
{
  // Instead of the previous ad-hoc method, we use Bareiss algorithm
  // to compute the determinant.
  // TODO - can do Cholesky, will be faster
  // !! TODO - Leibniz should be better when n <= 5 !?
  R sign = Math<R>::one();
  SquareMatrix<R, n+1> M;
  M(0,0) = Math<R>::one();
  // init
  for (size_t row = 0; row < n; row++)
    for (size_t col = 0; col < n; col++)
      M(row+1, col+1) = this->mat[row][col];
  for (size_t k = 1; k < n; k++) {
    if (M(k,k) == Math<R>::zero()) {
      bool found = false;
      for (size_t i = k+1; i <= n; i++) {
	if (M(i,k) != Math<R>::zero()) {
	  M.swap_rows(k,i);
	  sign = -sign;
	  found = true;
	  break;
	}
      }
      if (!found)
	return Math<R>::zero();
    }
    for (size_t i = k+1; i <= n; i++)
      for (size_t j = k+1; j <= n; j++)
	M(i,j) = (M(i,j)*M(k,k) - M(i,k)*M(k,j))/M(k-1,k-1);
  }
  return sign*M(n,n);
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
  SquareMatrix<R,n> inv;
  inv.set_identity();
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
  SquareMatrix<R,n> inv;
  inv.set_identity();
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
  L.set_identity();
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
  /*
#ifdef DEBUG
  // verify that L*Q*L^t = D
  SquareMatrix<R,n> diag;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      diag(i,j) = (i == j) ? D[i] : Math<R>::zero();
  assert(L*(*this)*L.transpose() == diag);
#endif
  */
  return true;
}

// !! TODO - This is integral Cholesky - should be able to
// determine which to call by the typename R
template<typename R, size_t n>
bool
SquareMatrix<R, n>::ldl(SquareMatrix<R, n>& L,  Vector<R,n> & D) const
{
  assert(is_symmetric());
  R prod_diag = Math<R>::one();
  R d, inner_sum;
  // This works but inefficiently - for some reason we get O(n^4) operations.
  // !! TODO - check it out later
  // Oh I see - we should do the L update in two passes...
  for (size_t i = 0; i < n; i++)
    {
      L(i, i) = prod_diag;
      d = prod_diag;
      for (size_t j = 0; j < i; j++)
	{
	  L(i, j) = Math<R>::zero();
	  for (size_t k = j; k < i; k++)
	    {
	      inner_sum = Math<R>::zero();
	      for (size_t r = 0; r <= k; r++)
		inner_sum += L(k, r)*((*this)(i,r))*L(k,j);
	      inner_sum *= -L(i, i) / D[k];
	      L(i,j) += inner_sum;
	    }
	  d = Math<R>::gcd(d, L(i, j));
	}
      for (size_t j = 0; j <= i; j++)
	L(i,j) /= d;
      D[i] = Math<R>::zero();
      for (size_t j = 0; j <= i; j++)
	for (size_t k = 0; k <= i; k++)
	  D[i] += L(i, j)*((*this)(j,k))*L(i, k);
      if (D[i] == Math<R>::zero()) return false;
      prod_diag = Math<R>::lcm(prod_diag, D[i]);
      for (size_t j = i+1; j < n; j++)
	L(i, j) = Math<R>::zero();
    }
  
#ifdef DEBUG
  // verify that L*Q*L^t = D
  SquareMatrix<R,n> diag;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      diag(i,j) = (i == j) ? D[i] : Math<R>::zero();
  assert(L*(*this)*L.transpose() == diag);
#endif
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
#ifdef DEBUG
  assert((row1 < n) && (row2 < n));
#endif 
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
#ifdef DEBUG
  assert((col1 < n) && (col2 < n));
#endif 
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
#ifdef DEBUG
  assert(row < n);
#endif  
  for (size_t col = 0; col < n; col++)
    mat[row][col] *= val;
  return;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::multiply_col(size_t col, const R & val)
{
#ifdef DEBUG
  assert(col < n);
#endif 
  for (size_t row = 0; row < n; row++)
    mat[row][col] *= val;
  return;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::add_row(size_t row_to, size_t row_from, const R & val)
{
#ifdef DEBUG
  assert((row_to < n) && (row_from < n));
#endif 
  for (size_t col = 0; col < n; col++) {
    mat[row_to][col] += val * mat[row_from][col];
  }
  return;
}

template<typename R, size_t n>
void SquareMatrix<R, n>::add_col(size_t col_to, size_t col_from, const R & val)
{
#ifdef DEBUG
  assert((col_to < n) && (col_from < n));
#endif 
  for (size_t row = 0; row < n; row++) {
    mat[row][col_to] += val * mat[row][col_from];
  }
  return;
}
  
// a general one, just in case
template<typename R, size_t n>
SquareMatrix<R, n> SquareMatrix<R, n>::inverse(void) const
{
#ifdef DEBUG
  assert(this->determinant() != Math<R>::zero());
#endif
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
  SquareMatrix<R, n> inv;
  inv.set_identity();
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
    if (max_val == Math<R>::zero()) {
      pivot_col++;
    }
    else {
      echelon.swap_rows(pivot_row, row_max);
      inv.swap_rows(pivot_row, row_max);
#ifdef DEBUG
      assert(inv*(*this) == echelon);
#endif
      R scalar = Math<R>::one() / echelon(pivot_row, pivot_col);
      echelon.multiply_row(pivot_row, scalar);
      inv.multiply_row(pivot_row, scalar);
#ifdef DEBUG
      assert(inv*(*this) == echelon);
#endif
      // for reduced row echelon form we need also the rows before
      for (size_t row = 0; row < pivot_row; row++) {
	factor = echelon(row, pivot_col);
	echelon(row, pivot_col) = Math<R>::zero();
	for (size_t col = pivot_col + 1; col < n; col++) {
	  echelon(row,col) -= factor * echelon(pivot_row, col);
	}
	for (size_t col = 0; col < n; col++) {
	  inv(row, col) -= factor * inv(pivot_row, col);
	}
#ifdef DEBUG
      assert(inv*(*this) == echelon);
#endif
      }
      for (size_t row = pivot_row+1; row < n; row++) {
	// factor = echelon(row,pivot_col) / echelon(pivot_row, pivot_col);
	factor = echelon(row, pivot_col);
	echelon(row, pivot_col) = Math<R>::zero();
	for (size_t col = pivot_col + 1; col < n; col++) {
	  echelon(row,col) -= factor * echelon(pivot_row, col);
	}
	for (size_t col = 0; col < n; col++) {
	  inv(row, col) -= factor * inv(pivot_row, col);
	}
#ifdef DEBUG
      assert(inv*(*this) == echelon);
#endif
      }
      
      pivot_row++;
      pivot_col++;
    }
  }
#ifdef DEBUG
  assert(inv*(*this) == echelon);
  assert((*this)*inv == identity());
#endif
  return inv;
}

template<typename R, size_t n>
SquareMatrix<R, n> SquareMatrix<R, n>::adjugate(size_t dim) const
{
  SquareMatrix<R, n> adj;
#ifdef DEBUG
  adj = SquareMatrix<R,n>::identity();
#endif
  // We will use this only for dim <= 4
  // and we write it down explicitly for each case
  // !! TODO - This is not the most effective way to do this
  const SquareMatrix<R,n> & a = (*this);
  if (dim <= n) {
    if (dim == 1) {
      adj(0,0) = Math<R>::one();
    }
    if (dim == 2) {
      adj(0,0) = a(1,1);
      adj(1,1) = a(0,0);
      adj(0,1) = -a(0,1);
      adj(1,0) = -a(1,0);
    }
    if (dim == 3) {
      adj(0,0) = a(1,1)*a(2,2) - a(1,2)*a(2,1);
      adj(1,0) = - a(1,0)*a(2,2) + a(1,2)*a(2,0);
      adj(2,0) = a(1,0)*a(2,1) - a(1,1)*a(2,0);
      adj(0,1) = - a(0,1)*a(2,2) + a(2,1)*a(0,2);
      adj(1,1) = a(0,0)*a(2,2) - a(0,2)*a(2,0);
      adj(2,1) = - a(0,0)*a(2,1) + a(0,1)*a(2,0);
      adj(0,2) = a(0,1)*a(1,2) - a(1,1)*a(0,2);
      adj(1,2) = - a(0,0)*a(1,2) + a(1,0)*a(0,2);
      adj(2,2) = a(1,1)*a(0,0) - a(1,0)*a(0,1);
    }
    if (dim == 4) {
      adj(0,0) = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2));
      adj(0,0) -= a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1));
      adj(0,0) += a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1));
      adj(0,1) = -a(1,0)*(a(2,2)*a(3,3)-a(2,3)*a(3,2));
      adj(0,1) += a(1,2)*(a(2,0)*a(3,3)-a(2,3)*a(3,0));
      adj(0,1) -= a(1,3)*(a(2,0)*a(3,2)-a(2,2)*a(3,0));
      adj(0,2) = a(1,0)*(a(2,1)*a(3,3)-a(2,3)*a(3,1));
      adj(0,2) -= a(1,1)*(a(2,0)*a(3,3)-a(2,3)*a(3,0));
      adj(0,2) += a(1,3)*(a(2,0)*a(3,1)-a(2,1)*a(3,0));
      adj(0,3) = -a(1,0)*(a(2,1)*a(3,2)-a(2,2)*a(3,1));
      adj(0,3) += a(1,1)*(a(2,0)*a(3,2)-a(2,2)*a(3,0));
      adj(0,3) -= a(1,2)*(a(2,0)*a(3,1)-a(2,1)*a(3,0));
      adj(1,0) = -a(0,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2));
      adj(1,0) += a(2,1)*(a(0,2)*a(3,3)-a(3,2)*a(0,3));
      adj(1,0) -= a(3,1)*(a(0,2)*a(2,3)-a(2,2)*a(0,3));
      adj(1,1) = a(0,0)*(a(2,2)*a(3,3)-a(2,3)*a(3,2));
      adj(1,1) -= a(0,2)*(a(2,0)*a(3,3)-a(2,3)*a(3,0));
      adj(1,1) += a(0,3)*(a(2,0)*a(3,2)-a(2,2)*a(3,0));
      adj(1,2) = -a(0,0)*(a(2,1)*a(3,3)-a(2,3)*a(3,1));
      adj(1,2) += a(0,1)*(a(2,0)*a(3,3)-a(3,0)*a(2,3));
      adj(1,2) -= a(0,3)*(a(2,0)*a(3,1)-a(3,0)*a(2,1));
      adj(1,3) = a(0,0)*(a(2,1)*a(3,2)-a(3,1)*a(2,2));
      adj(1,3) -= a(0,1)*(a(2,0)*a(3,2)-a(3,0)*a(2,2));
      adj(1,3) += a(0,2)*(a(2,0)*a(3,1)-a(2,1)*a(3,0));
      adj(2,0) = a(0,1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3));
      adj(2,0) -= a(1,1)*(a(0,2)*a(3,3)-a(3,2)*a(0,3));
      adj(2,0) += a(3,1)*(a(0,2)*a(1,3)-a(1,2)*a(0,3));
      adj(2,1) = -a(0,0)*(a(1,2)*a(3,3)-a(3,2)*a(1,3));
      adj(2,1) += a(1,0)*(a(0,2)*a(3,3)-a(0,3)*a(3,2));
      adj(2,1) -= a(3,0)*(a(0,2)*a(1,3)-a(0,3)*a(1,2));
      adj(2,2) = a(0,0)*(a(1,1)*a(3,3)-a(1,3)*a(3,1));
      adj(2,2) -= a(0,1)*(a(1,0)*a(3,3)-a(1,3)*a(3,0));
      adj(2,2) += a(0,3)*(a(1,0)*a(3,1)-a(1,1)*a(3,0));
      adj(2,3) = -a(0,0)*(a(1,1)*a(3,2)-a(1,2)*a(3,1));
      adj(2,3) += a(0,1)*(a(1,0)*a(3,2)-a(3,0)*a(1,2));
      adj(2,3) -= a(0,2)*(a(1,0)*a(3,1)-a(3,0)*a(1,1));
      adj(3,0) = -a(0,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3));
      adj(3,0) += a(1,1)*(a(0,2)*a(2,3)-a(2,2)*a(0,3));
      adj(3,0) -= a(2,1)*(a(0,2)*a(1,3)-a(1,2)*a(0,3));
      adj(3,1) = a(0,0)*(a(1,2)*a(2,3)-a(1,3)*a(2,2));
      adj(3,1) -= a(1,0)*(a(0,2)*a(2,3)-a(0,3)*a(2,2));
      adj(3,1) += a(2,0)*(a(0,2)*a(1,3)-a(1,2)*a(0,3));
      adj(3,2) = -a(0,0)*(a(1,1)*a(2,3)-a(2,1)*a(1,3));
      adj(3,2) += a(1,0)*(a(0,1)*a(2,3)-a(0,3)*a(2,1));
      adj(3,2) -= a(2,0)*(a(0,1)*a(1,3)-a(0,3)*a(1,1));
      adj(3,3) = a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1));
      adj(3,3) -= a(0,1)*(a(1,0)*a(2,2)-a(1,2)*a(2,0));
      adj(3,3) += a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));
    }
  }
#ifdef DEBUG
  SquareMatrix<R, n> a_copy = SquareMatrix<R,n>::identity();
  for (size_t row = 0; row < dim; row++)
    for (size_t col = 0; col < dim; col++)
      a_copy(row,col) = a(row,col);
  R det = a_copy.determinant();
  SquareMatrix<R,n>  prod = adj*a;
  for (size_t row = 0; row < dim; row++)
    for (size_t col = 0; col < dim; col++)
      assert(prod(row,col) == ((row==col) ? det : Math<R>::zero()));
#endif
  return adj;
}

// static functions
template<typename R, size_t n>
Rational<R>
SquareMatrix<R,n>::inner_product(const SquareMatrix<R, n> & F,
				 const SquareMatrix<Rational<R>, n> & S,
				 size_t idx1, size_t idx2)
{
#ifdef DEBUG
  assert((idx1 < n) && (idx2 < n));
#endif 
  Rational<R> ans = Math< Rational<R> >::zero();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ans += S(idx1, i) * F(i,j) * S(idx2, j);
  return ans;
}

// Here we compute a highly specialized hermite form
// Assuming the matrix is non-singular
// and that we are looking for the hermite form of the matrix
// concatenated to the scalar matrix d*I
template<typename R, size_t n>
SquareMatrix<R,n> SquareMatrix<R, n>::hermite_form(const R & d) const
{
  R a,b,x,y,g,q_a,q_b;
  SquareMatrix<R, n> H = d*SquareMatrix<R,n>::identity();
  for (size_t row = 0; row < n; row++) {
    Vector<R, n> b_primes = (*this)[row];
    // Here we compute the HNF of H and this row and store it in H
    for (size_t pivot = 0; pivot < n; pivot++) {
      a = H(pivot,pivot);
      b = b_primes[pivot];
      g = Math<R>::xgcd(a,b,x,y);
      q_a = a / g;
      q_b = b / g;
      Vector<R, n> g_h_prime = x*H[pivot] + y*b_primes;
      b_primes = q_a*b_primes-q_b*H[pivot];
      for (size_t j = pivot; j < n; j++) {
	R scalar = b_primes[j] / H(j,j);
	b_primes -= scalar*H[j];
      }
      for (size_t col = 0; col < n; col++)
	H(pivot, col) = g_h_prime[col];
    }
    for (size_t pivot = n-1; pivot > 0; pivot--) {
      for (size_t col = pivot; col < n; col++) { 
	R q = H(pivot-1,col) / H(col, col);
	for (size_t j = col; j < n; j++)
	  H(pivot-1, j) -= q*H(col, j);
      }
    }
  }
  return H;
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

template<typename R, size_t n>
void SquareMatrix<R,n>::set_identity(void)
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      (*this)(i,j) = (i == j) ? Math<R>::one() : Math<R>::zero();
  
  return;
}

template<typename R, typename S, size_t n>
void SquareMatrixFp<R, S, n>::set_identity(void)
{
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      (*this)(i,j) = (i == j) ? Math<R>::one() : Math<R>::zero();
      (*this)(i,j).set_field(this->GF);
    }
  
  return;
}

template<typename R, typename S, size_t n>
SquareMatrixFp<R, S, n>
SquareMatrixFp<R, S, n>::operator*(const SquareMatrixFp<R, S, n>& other) const
{
  SquareMatrixFp<R, S, n> prod(this->GF);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      prod(i,j) = 0;
      for (size_t k = 0; k < n; k++)
	prod(i,j) += this->mat[i][k]*other(k,j);
    }
  return prod;
}
