// Implementation of templated functions from QuadForm.h

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
