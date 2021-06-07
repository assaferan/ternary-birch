template<typename R>
NumberField<R>::NumberField(const UnivariatePoly<R> & mod)
{
  for (int i = 0; i <= mod.degree(); i++)
    this->f += mod.coefficient(i)*UnivariatePoly< Rational<R> >::x(i);
  
}

template<typename R>
NumberFieldElement<R> & NumberFieldElement<R>::operator=(const R & a)
{
  this->elt = a;
  return (*this);
}

template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator-() const
{
  NumberFieldElement<R> neg(this->K);
  neg.elt = -(this->elt);
  
  return neg;
}

template<typename R>
NumberFieldElement<R>
NumberFieldElement<R>::operator+(const NumberFieldElement<R> & other) const
{
  NumberFieldElement<R> sum(this->K);
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
  NumberFieldElement<R> prod(this->K);
  prod.elt = (this->elt)  * other.elt % (this->K->modulus());

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
  NumberFieldElement<R> prod(this->K);
  prod.elt = (this->elt) * a;

  return prod;
}

template<typename R>
NumberFieldElement<R> NumberFieldElement<R>::operator/(const R & a) const
{
  NumberFieldElement<R> quo(this->K);
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
  UnivariatePoly< Rational<R> >::xgcd(this->elt, this->K->modulus(), s, t);
  
  NumberFieldElement<R> inv(this->K, s);
  return inv;
}






