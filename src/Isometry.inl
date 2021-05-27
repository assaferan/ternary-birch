// implementation for Isometry.h

// TODO - change everythign from QuadForm to SquareMatrix

template<typename R, size_t n>
SquareMatrix<R, n>
Isometry<R,n>::transform(const SquareMatrix<R, n>& from) const
{
  return (this->a).transpose()*from*(this->a) / (this->scale * this->scale);
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
	  val += this->a(k,i)*from.bilinear_form()(k,l)*
	    this->a(l,j);
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

template<typename R, size_t n>
Isometry<R, n> Isometry<R,n>::inverse(void) const
{
  // !! TODO - should be able to invert without using rationals
  // for example, can always track back (save the inverse for the ride)
  SquareMatrix< Rational<R>, n> a_rat;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      a_rat(i,j) = this->a(i,j);
      a_rat(i,j) /= this->scale;
    }
  a_rat = a_rat.inverse();
  SquareMatrix<R, n> a_inv;
  // Since this is an isometry, the inverse should be integral
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      a_inv(i,j) = (this->scale * a_rat(i,j)).floor();
#ifdef DEBUG
  R scale2 = (this->scale)*(this->scale);
  assert((a_inv * (this->a) == scale2*SquareMatrix<R, n>::identity()));
#endif
  return Isometry(a_inv, this->scale);
}

template<typename R, size_t n>
void Isometry<R,n>::rescale(void)
{
  R d = this->scale;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      d = Math<R>::gcd(d, this->a(i,j));
  
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->a(i,j) /= d;

  this->scale /= d;
  return;
}
