// For implementations

template<typename R>
template <size_t n>
Matrix<R>::Matrix(const R data[n][n])
  : nrows_(n), ncols_(n), data_(n*n)
{ size_t idx = 0;
  for (size_t row = 0; row < nrows_; row++)
    for (size_t col = 0; col < ncols_; col++)
      data_[idx++] = data[row][col];
}

template<typename R>
template <size_t n>
Matrix<R>::Matrix(const SquareMatrix<R, n> & mat)
  : nrows_(n), ncols_(n), data_(n*n)
{
  size_t idx = 0;
  for (size_t row = 0; row < nrows_; row++)
    for (size_t col = 0; col < ncols_; col++)
      data_[idx++] = mat(row, col);
}

template<typename R>
R Matrix<R>::determinant() const
{		
  assert(nrows_ == ncols_);
  size_t n = nrows_;
  Matrix<R> M(n+1, n+1);
  M(0,0) = Math<R>::one();
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

template<typename R>
void Matrix<R>::swap_rows(size_t row1, size_t row2)
{
  R tmp;
  for (size_t col = 0; col < this->ncols_; col++) {
    tmp = (*this)(row1, col);
    (*this)(row1, col) = (*this)(row2, col);
    (*this)(row2, col) = tmp;
  }
  return;
}

// in place echelon form, returns the rank and trans is the transformation
template<typename R>
size_t Matrix<R>::row_echelon(Matrix<R> & echelon, Matrix<R>& trans)
{
  trans = identity(echelon.nrows());
  size_t pivot_row = 0;
  size_t pivot_col = 0;
  size_t row_max;
  int max_val, val;
  R factor;
  
  while ((pivot_row < echelon.nrows()) && (pivot_col < echelon.ncols())) {
    row_max = pivot_row;
    max_val = abs(echelon(row_max, pivot_col));
    for (size_t row = pivot_row+1; row < echelon.nrows(); row++) {
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
      trans.swap_rows(pivot_row, row_max);
      for (size_t row = pivot_row+1; row < echelon.nrows(); row++) {
	factor = echelon(row,pivot_col) / echelon(pivot_row, pivot_col);
	echelon(row, pivot_col) = 0;
	for (size_t col = pivot_col + 1; col < echelon.ncols(); col++) {
	  echelon(row,col) -= factor * echelon(pivot_row, col);
	}
	for (size_t col = 0; col < echelon.ncols(); col++) {
	  trans(row,col) -= factor * trans(pivot_row, col);
	}
      }
      
      pivot_row++;
      pivot_col++;
    }
  }
  return pivot_col;
}

template<typename R>
size_t Matrix<R>::rank() const
{  
  Matrix<R> echelon((*this));
  Matrix<R> trans(echelon.nrows(), echelon.nrows());
  return row_echelon(echelon, trans);
}

template<typename R, typename S>
size_t MatrixFp<R, S>::rank() const
{  
  MatrixFp<R, S> echelon((*this));
  MatrixFp<R, S> trans(this->GF, echelon.nrows(), echelon.nrows());
  return Matrix<FpElement<R,S> >::row_echelon(echelon, trans);
}

template<typename R>
Matrix<R> Matrix<R>::kernel() const {
  return this->transpose().left_kernel();
}

template<typename R>
Matrix<R> Matrix<R>::left_kernel() const {
  Matrix<R> echelon((*this));
  Matrix<R> trans(echelon.nrows(), echelon.nrows());
  size_t rank = row_echelon(echelon, trans);
  // getting the zero rows
  Matrix<R> kernel(nrows_ - rank, ncols_);
   for (size_t row = rank; row < nrows_; row++)
    for (size_t col = 0; col < ncols_; col++)
      kernel(row-rank,col) = trans(row, col);
  return kernel;
}

template<typename R>
Matrix<R> Matrix<R>::diagonal_join(const std::vector< Matrix<R> > & mats)
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

template<typename R>
Matrix<R> Matrix<R>::identity(size_t n) {
  Matrix<R> id(n,n);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      id(i,j) = (i == j) ? Math<R>::one() : Math<R>::zero();
  return id;
}

// TODO - just change access resolution to the same vector instead
template<typename R>
Matrix<R> Matrix<R>::transpose() const
{
  std::vector<R> data_t(nrows_*ncols_);
  size_t idx = 0;
  for (size_t col = 0; col < ncols_; col++)
    for (size_t row = 0; row < nrows_; row++)
      data_t[idx++] = (*this)(row, col);
  Matrix<R> tr(data_t, ncols_, nrows_);
  return tr;
}

template<typename R>
Matrix<R> Matrix<R>::operator*(const Matrix<R> & other) const
{
  size_t nrows = this->nrows_;
  size_t ncols = other.ncols_;
  assert( this->ncols_ == other.nrows_ );
  Matrix<R> prod(nrows, ncols);
    
  for (size_t row = 0; row < nrows; row++)
    for (size_t col = 0; col < ncols; col++) {
      prod(row,col) = Math<R>::zero();
      for (size_t j = 0; j < this->ncols_; j++)
	prod(row,col) += (*this)(row,j)*other(j,col);
    }
  return prod;
}
