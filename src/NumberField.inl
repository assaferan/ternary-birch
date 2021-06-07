template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator-() const
{
  NumberFieldElement<R> neg;
  neg.elt = -(this->elt);
  
  return neg;
}

template<typename R>
NumberFieldElement<R>
NumberFieldElement<R>::operator+(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> sum;
  sum.elt = (this->elt) + other.elt;

  return sum;
}

template<typename R>
NumberFieldElement<R>
NumberFieldElement<R>::operator-(const NumberFieldElement<R> & other) const
{
  return (*this)+(-other);
}

template<typename R>
NumberFieldElement<R>
NumberFieldElement<R>::operator*(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> prod;
  prod.elt = (this->elt)  * other.elt % NumberFieldElement<R>::modulus;

  return prod;
}

template<typename R>
NumberFieldElement<R>
NumberFieldElement<R>::operator/(const NumberFieldElement<R> & other) const
{
  return (*this) * other.inverse();
}

template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator*(const R & a) const
{
  NumberFieldElement<R> prod;
  prod.elt = (this->elt) * a;

  return prod;
}

template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator/(const R & a) const
{
  NumberFieldElement<R> quo;
  quo.elt = (this->elt) / a;

  return quo;
}

template<typename R>
NumberFieldElement<R> &
NumberFieldElement<R>::operator+=(const NumberFieldElement<R> & other)
{
  this->elt += other.elt;
  return (*this);
}

template<typename R>
NumberFieldElement<R> &
NumberFieldElement<R>::operator-=(const NumberFieldElement<R> & other)
{
  this->elt -= other.elt;
  return (*this);
}

template<typename R>
NumberFieldElement<R> &
NumberFieldElement<R>::operator*=(const NumberFieldElement<R> & other)
{
  this = (*this)*other;
  return (*this);
}

template<typename R>
NumberFieldElement<R> &
NumberFieldElement<R>::operator/=(const NumberFieldElement<R> & other)
{
  this = (*this)/other;
  return (*this);
}

template<typename R>
NumberFieldElement<R>& NumberFieldElement<R>::operator*=(const R & a)
{
  this->elt *= a;
  return (*this);
}

template<typename R>
NumberFieldElement<R>& NumberFieldElement<R>::operator/=(const R & a)
{
  this->elt /= a;
  return (*this);
}

template<typename R>
NumberFieldElement<R>
NumberFieldElement<R>::inverse() const
{
  UnivariatePoly< Rational<R> > s,t;
  xgcd(this->elt, NumberFieldElement<R>::modulus, s, t);
  
  NumberFieldElement<R> inv(s);
  return inv;
}






