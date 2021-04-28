template <typename R>
Rational<R> Rational<R>::operator+(const Rational<R> &other) const
{
  Rational<R> sum;
  sum.num_ = (this->num_) * other.denom_ + (this->denom_) * other.num_;
  sum.denom_ = (this->denom_) * other.denom_;
  sum.reduce();
  return sum;
}

template <typename R>
Rational<R> Rational<R>::operator*(const Rational<R> &other) const
{
  Rational<R> prod;
  prod.num_ = (this->num_) * other.num_;
  prod.denom_ = (this->denom_) * other.denom_;
  prod.reduce();
  return prod;
}

template <typename R>
Rational<R> Rational<R>::operator/(const Rational<R> &other) const
  {
  Rational<R> prod;
  prod.num_ = (this->num_) * other.denom__;
  prod.denom_ = (this->denom_) * other.num_;
  prod.reduce();
  return prod;
}

template <typename R>
bool Rational<R>::operator==(const Rational<R> &other) const
{
  return (this->num_) * other.denom_ == other.num_ * (this->denom_); 
}

template <typename R>
bool Rational<R>::operator<(const Rational<R> &other) const
{
  return (this->num_) * other.denom_ < other.num_ * (this->denom_); 
}

// can also implement a general gcd - maybe later
template class Rational<Z>;

void Rational<Z>::reduce(void)
{
  Z d;
  mpz_gcd(&d, num_, denom_);
  num_ /= d;
  denom_ /= d;
  return;
}