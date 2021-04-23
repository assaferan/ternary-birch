// Implementation of templated functions from QuadForm.h
#include <type_traits>

// c-tors
template<typename R, size_t Rank>
QuadForm<R, Rank>::QuadForm(const RVec& coeffs)
{
  size_t idx = 0;
  for (size_t row = 0; row < Rank; row++)
    {
      for (size_t col = 0; col < row; col++)
	{	
	  this->B_[row][col] = coeffs[idx];
	  this->B_[col][row] = coeffs[idx++];
	}
      this->B_[row][row] = coeffs[idx++];
    }
}

template<typename R, size_t Rank>
R QuadForm<R, Rank>::discriminant(void) const
{
  if (Rank == 3) 
        return this->a_ * (4 * this->b_ * this->c_ - this->f_ * this->f_) -
            this->b_ * this->g_ * this->g_ +
            this->h_ * (this->f_ * this->g_ - this->c_ * this->h_);
  else
  {
  // Instead of the previous ad-hoc method, we use Bareiss algorithm
  // to compute the determinant.
  // TODO - can do Cholesky, will be faster
  // !! TODO - Leibniz should be better when Rank <= 5
    R M[Rank+1][Rank+1];
    M[0][0] = 1;
    // init
    for (size_t row = 0; row < Rank; row++)
      for (size_t col = 0; col < Rank; col++)
        M[row+1][col+1] = this->B_[row][col];
    for (size_t k = 1; k < Rank; k++)
      for (size_t i = k+1; i <= Rank; i++)
        for (size_t j = k+1; j <= Rank; j++)
          M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j])/M[k-1][k-1];
    if (Rank % 2 == 0)	  
      return M[Rank][Rank];
    else
      return M[Rank][Rank]/2;
  }
}
