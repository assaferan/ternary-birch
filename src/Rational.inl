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
  prod.num_ = (this->num_) * other.denom_;
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

template<typename R>
R my_gcd(const R & a, const R & b)
{
  if (b < 0) return my_gcd(a,-b);
  if (a < 0) return my_gcd(-a,b);
  if (a < b) return my_gcd(b,a);
  if (b == 0) return a;
  return my_gcd(b, a % b);
}

template<typename R>
void Rational<R>::reduce(void)
{
  R d = my_gcd(num_, denom_);
  num_ /= d;
  denom_ /= d;
  return;
}
