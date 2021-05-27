// implementation for Isometry.h

// TODO - change everythign from QuadForm to SquareMatrix

template<typename R, size_t n>
SquareMatrix<R, n>
Isometry<R,n>::transform(const SquareMatrix<R, n>& from, R scalar) const
{
  // Is this the right direction? check carefully
  return (this->a).transpose()*from*(this->a) / scalar;
}

template<typename R, size_t n>
bool Isometry<R,n>::is_isometry(const QuadForm<R, n>& from,
				const QuadForm<R, n>& to,
				R scalar) const
{
  R val;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      val = 0;
      for (size_t k = 0; k < n; k++)
	for (size_t l = 0; l < n; l++)
	  val += this->a(i,k)*from.bilinear_form()(k,l)*
	    this->a(j,l);
      if (val != to.bilinear_form()(i,j) * scalar)
	return false;
    }
  return true;
}

template<typename R, size_t n>
void Isometry<R,n>::update_perm(const Vector<size_t, n> & perm) {
  SquareMatrix<R,n> temp;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      temp(i,j) = this->a(i,j);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->a(perm[i],j) = temp(i,j);
  return;
}
