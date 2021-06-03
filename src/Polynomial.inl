// create the constant polynomial
template<typename R>
UnivariatePoly<R>::UnivariatePoly(const R & a)
{
  if (a != Math<R>::zero()) {
    this->coeffs.resize(1);
    this->coeffs[0] = a;
  }
}

// create polynomial from coefficients
template<typename R>
UnivariatePoly<R>::UnivariatePoly(const std::vector<R> & vec)
  : coeffs(vec)
{}

// create the polynomial x
template<typename R>
UnivariatePoly<R> UnivariatePoly<R>::x()
{
  UnivariatePoly<R> p;
  p.coeffs.resize(2);
  p[0] = Math<R>::zero();
  p[1] = Math<R>::one();
  return p;
}

// coefficient of x^i
template<typename R>
R UnivariatePoly<R>::coefficient(size_t i) const
{
  if (i < this->coeffs.size())
    return this->coeffs[i];
  return Math<R>::zero();
}

// conversion, assignment operator
template<typename R>
UnivariatePoly<R> &
UnivariatePoly<R>::operator=(const UnivariatePoly<R> & other)
{
  if (this != (&other)) {
    this->coeffs = other.coeffs;
  }
  return (*this); 
}

template<typename R>
UnivariatePoly<R> & UnivariatePoly<R>::operator=(const R & a)
{
  this->coeffs.resize(1);
  this->coeffs[0] = a;
  return (*this);
}
  
// arithmetic
template<typename R>
UnivariatePoly<R> UnivariatePoly<R>::operator-() const
{
  UnivariatePoly<R> neg;
  neg.coeffs.resize(this->coeffs.size());
  for (size_t i = 0; i < this->coeffs.size(); i++)
    neg.coeffs[i] = -this->coeffs[i];
  return neg;
}

template<typename R>
UnivariatePoly<R>
UnivariatePoly<R>::operator+(const UnivariatePoly<R> & other) const
{
  UnivariatePoly<R> sum;
  if (this->coeffs.size() < other.coeffs.size())
    return other + (*this);
  // here we may assume this is the polynomial with the larger degree
  sum.coeffs.resize(this->coeffs.size());
  size_t i;
  for (i = 0; i < other.coeffs.size(); i++)
    sum.coeffs[i] = this->coeffs[i] + other.coeffs[i];
  for (; i < this->coeffs.size(); i++)
    sum.coeffs[i] = this->coeffs[i];

  // eliminate redundant zeros
  if (this->coeffs.size() == other.coeffs.size()) {
    i = sum.coeffs.size();
    while((i > 0) && (sum.coeffs[i-1] == Math<R>::zero())) i--;
    sum.coeffs.resize(i);
  }
  return sum;
}

template<typename R>
UnivariatePoly<R>
UnivariatePoly<R>::operator-(const UnivariatePoly<R> & other) const
{
  return (*this) + (-other);
}

template<typename R>
UnivariatePoly<R>
UnivariatePoly<R>::operator*(const UnivariatePoly<R> & other) const
{
  UnivariatePoly<R> prod;
  prod.coeffs.resize(this->degree()+other.degree()+1);
  std::fill(prod.coeffs.begin(), prod.coeffs.end(), Math<R>::zero());
  size_t i, j;
  for (i = 0; i < this->coeffs.size(); i++)
    for (j = 0; j < other.coeffs.size(); j++)
      prod.coeffs[i+j] += this->coeffs[i] * other.coeffs[j];
  
  return prod;
}

template<typename R>
UnivariatePoly<R> UnivariatePoly<R>::operator*(const R & a) const
{
  UnivariatePoly<R> prod;
  prod.coeffs.resize(this->coeffs.size());
  for (size_t i = 0; i < this->coeffs.size(); i++)
    prod.coeffs[i] = a * this->coeffs[i];
  
  return prod;
}

template<typename R>
UnivariatePoly<R> &
UnivariatePoly<R>::operator+=(const UnivariatePoly<R> & other)
{
  size_t i, deg;
  deg = this->coeffs.size();
  if (deg < other.coeffs.size()) {
    this->coeffs.resize(other.coeffs.size());
    for (i = 0; i < deg; i++)
      this->coeffs[i] += other.coeffs[i];
    for (; i < other.coeffs.size(); i++)
      this->coeffs[i] = other.coeffs[i];
  }
  else {
    for (i = 0; i < other.coeffs.size(); i++)
      this->coeffs[i] += other.coeffs[i];
     // eliminate redundant zeros
    if (this->coeffs.size() == other.coeffs.size()) {
      while((deg > 0) && (this->coeffs[deg-1] == Math<R>::zero())) deg--;
      this->coeffs.resize(deg);
    }
  }
  return (*this);
}

template<typename R>
UnivariatePoly<R> &
UnivariatePoly<R>::operator-=(const UnivariatePoly<R> & other)
{
  return ((*this) += (-other));
}

template<typename R>
UnivariatePoly<R> &
UnivariatePoly<R>::operator*=(const UnivariatePoly<R> & other)
{
  (*this) = (*this)*other;
  return (*this);
}

template<typename R>
UnivariatePoly<R>& UnivariatePoly<R>::operator*=(const R & a)
{
  for (size_t i = 0; i < this->coeffs.size(); i++)
    this->coeffs[i] *= a;
  
  return (*this);
}

template<typename R>
UnivariatePoly<R>
UnivariatePoly<R>::evaluate(const UnivariatePoly<R> & f) const
{
  UnivariatePoly<R> comp;
  UnivariatePoly<R> f_i(Math<R>::one());
  for (size_t i = 0; i < this->coeffs.size(); i++) {
    comp += this->coeffs[i]*f_i;
    f_i *= f;
  }
  return comp;
}

template<typename R>
R UnivariatePoly<R>::evaluate(const R & a) const
{
  R res;
  R a_i = Math<R>::one();
  for (size_t i = 0; i < this->coeffs.size(); i++) {
    res += this->coeffs[i]*a_i;
    a_i *= a;
  }
  return res;
}

// booleans
template<typename R>
bool UnivariatePoly<R>::operator==(const UnivariatePoly<R> & other) const
{
  if (this->coeffs.size() != other.coeffs.size())
    return false;
  for (size_t i = 0; i < this->coeffs.size(); i++)
    if(this->coeffs[i] != other.coeffs[i])
      return false;
  
  return true;
}

template<typename R>
bool UnivariatePoly<R>::operator!=(const UnivariatePoly<R> & other) const
{
  return !((*this)==other);
}

template<typename R>
bool UnivariatePoly<R>::operator==(const R & a) const
{
  if (a == Math<R>::zero()) return this->is_zero();
  if (this->coeffs.size() != 1)
    return false;

  return (this->coeffs[0] == a);
}

template<typename R>
bool UnivariatePoly<R>::operator!=(const R & a) const
{
  return !((*this)==a);
}

template<typename R>
std::ostream& operator<<(std::ostream& os, const UnivariatePoly<R> & p)
{
  size_t deg = p.degree();
  for (size_t i = deg+1; i > 0; i--) {
    R coeff = p.coefficient(i-1);
    if (coeff != Math<R>::zero()) {
      if ((i <= deg) && (coeff > Math<R>::zero()))
	os << "+";
      if (coeff != Math<R>::one())
	os << coeff;
      if (i > 1)
	os << "x";
      if (i > 2)
	os << "^" << (i-1);
    }
  }
  return os;
}

template<typename R>
std::vector< std::pair< UnivariatePoly<R>, size_t > >
UnivariatePolynomial<R>::factor() const
{
  std::vector< std::pair< UnivariatePoly<R>, size_t > > fac;
  
  return fac;
}

// PolynomialFp

// create the zero polynomial
template<typename R, typename S>
PolynomialFp<R,S>::PolynomialFp(std::shared_ptr<const Fp<R,S>> GF)
{
  this->GF = GF;
}

// create the constant polynomial 
template<typename R, typename S>
PolynomialFp<R,S>::PolynomialFp(const FpElement<R, S> & a)
{
  this->GF = a.field();
  std::multiset<size_t> empty_set;
  this->mons[empty_set] = a;
}

// create the constant polynomial 
template<typename R, typename S>
PolynomialFp<R,S>::PolynomialFp(std::shared_ptr<const Fp<R,S>> GF, 
				  const R & a)
{
  this->GF = GF;
  std::multiset<size_t> empty_set;
  FpElement<R,S> elt(GF, a);
  this->mons[empty_set] = elt;
}

// create the polynomial x_i
template<typename R, typename S>
PolynomialFp<R,S>
PolynomialFp<R,S>::x(std::shared_ptr<const Fp<R,S>> GF, size_t i)
{
  PolynomialFp<R,S> x_i(GF);

  std::multiset<size_t> singleton;
  singleton.insert(i);
  FpElement<R,S> one(GF,1);
  x_i.mons[singleton] = one;

  return x_i;
}

template<typename R, typename S>
template<size_t n>
PolynomialFp<R,S>::PolynomialFp(const SquareMatrixFp<R, S, n> & q)
{
  this->GF = q.field();
  
  for (size_t i = 0; i < n; i++)
    for (size_t j = i; j < n; j++) {
      std::multiset<size_t> mon;
      mon.insert(i);
      mon.insert(j);
      this->mons[mon] = q(i,j);
    }
  if (this->GF->prime() != 2) {
    FpElement<R,S> two(GF,2);
    for (size_t i = 0; i < n; i++) {
      std::multiset<size_t> mon;
      mon.insert(i);
      mon.insert(i);
      this->mons[mon] /= two;
    }
  }
}

// returns the coefficient of a monomial
template<typename R, typename S>
FpElement<R, S>
PolynomialFp<R,S>::coefficient(const std::multiset<size_t> & mon) const
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  
  it = this->mons.find(mon);
  if (it == mons.end()) {
    FpElement<R,S> zero(this->GF,0);
    return zero;
  }
  return it->second;
}

template<typename R, typename S>
FpElement<R, S> PolynomialFp<R,S>::const_coefficient() const {
  std::multiset<size_t> empty_set;
  return this->coefficient(empty_set);
}

// coefficient of x_i
template<typename R, typename S>
FpElement<R, S> PolynomialFp<R,S>::coefficient(size_t i) const {
  std::multiset<size_t> singleton;
  singleton.insert(i);
  return this->coefficient(singleton);
}

// coefficient of x_i x_j
template<typename R, typename S>
FpElement<R, S> PolynomialFp<R,S>::coefficient(size_t i, size_t j) const
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  std::multiset<size_t> mon;
  mon.insert(i);
  mon.insert(j);
  return this->coefficient(mon);
}

template<typename R, typename S>
PolynomialFp<R,S> PolynomialFp<R,S>::quadratic_part() const {
  PolynomialFp<R,S> quad(this->GF);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->mons.begin(); i != this->mons.end(); i++) {
    if ((i->first).size() == 2)
      quad.mons[i->first] = i->second;
  }
  
  return quad;
}

template<typename R, typename S>
std::vector< FpElement<R,S> > PolynomialFp<R,S>::linear_part(size_t rank) const
{
  std::vector< FpElement<R,S> > linear;
  for (size_t i = 0; i < rank; i++)
    linear.push_back(this->coefficient(i));
  return linear;
}

template<typename R, typename S>
int PolynomialFp<R,S>::degree(size_t i) const
{
  int deg = -1;
  int mon_deg;
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  
  for (it = this->mons.begin(); it != this->mons.end(); it++) {
    if (it->second != 0) {
      mon_deg = it->first.count(i);
      if (deg < mon_deg)
	deg = mon_deg;
    }
  }

  return deg;
}

template<typename R, typename S>
PolynomialFp<R,S> & PolynomialFp<R,S>::operator=(const PolynomialFp<R,S> & other)
{
  if (this != (&other)) {
    this->GF = other.GF;
    this->mons = other.mons;
  }
  return (*this);
}

template<typename R, typename S>
PolynomialFp<R,S> & PolynomialFp<R,S>::operator=(const FpElement<R,S> & a)
{
  this->GF = a.field();
  this->mons.clear();
  std::multiset<size_t> empty_set;
  this->mons[empty_set] = a;
  return (*this);
}

template<typename R, typename S>
PolynomialFp<R,S> PolynomialFp<R,S>::operator-() const
{
  PolynomialFp<R,S> neg(this->GF);
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  for (it = this->mons.begin(); it != this->mons.end(); it++) {
    neg.mons[it->first] = -it->second;
  }
  
  return neg;
}

template<typename R, typename S>
PolynomialFp<R,S> PolynomialFp<R,S>::operator+(const PolynomialFp<R,S> & other) const
{
  PolynomialFp<R,S> sum(this->GF);
  FpElement<R,S> zero(this->GF, 0);
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  for (it = this->mons.begin(); it != this->mons.end(); it++) {
    sum.mons[it->first] = it->second;
  }
  for (it = other.mons.begin(); it != other.mons.end(); it++) {
    it2 = sum.mons.find(it->first);
    if (it2 == sum.mons.end())
      sum.mons[it->first] = zero;
    sum.mons[it->first] += it->second;
  }
  
  return sum;
}

template<typename R, typename S>
PolynomialFp<R,S> PolynomialFp<R,S>::operator-(const PolynomialFp<R,S> & other) const
{
  PolynomialFp<R,S> diff(this->GF);
  FpElement<R,S> zero(this->GF, 0);
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  for (it = this->mons.begin(); it != this->mons.end(); it++) {
    diff.mons[it->first] = it->second;
  }
  for (it = other.mons.begin(); it != other.mons.end(); it++) {
    it2 = diff.mons.find(it->first);
    if (it2 == diff.mons.end())
      diff.mons[it->first] = zero;
    diff.mons[it->first] -= it->second;
  }
  
  return diff;
}

template<typename R, typename S>
PolynomialFp<R,S> PolynomialFp<R,S>::operator*(const FpElement<R,S> & a) const
{
  PolynomialFp<R,S> prod(this->GF);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  for (it = this->mons.begin(); it != this->mons.end(); it++) {
    prod.mons[it->first] = a*it->second;
  }
  
  return prod;
}

template<typename R, typename S>
PolynomialFp<R,S> PolynomialFp<R,S>::operator*(const PolynomialFp<R,S> & other) const
{
  PolynomialFp<R,S> prod(this->GF);
  FpElement<R,S> zero(this->GF, 0);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i,j,loc;
  for (i = this->mons.begin(); i != this->mons.end(); i++)
    for (j = other.mons.begin(); j != other.mons.end(); j++) {
      std::multiset<size_t> mon;
      mon.insert(i->first.begin(), i->first.end());
      mon.insert(j->first.begin(), j->first.end());
      loc = prod.mons.find(mon);
      if (loc == prod.mons.end())
	prod.mons[mon] = zero;
      prod.mons[mon] += (i->second)*(j->second);
  }
  
  return prod;
}

template<typename R, typename S>
PolynomialFp<R,S> & PolynomialFp<R,S>::operator+=(const PolynomialFp<R,S> & other)
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  
  for (it = other.mons.begin(); it != other.mons.end(); it++) {
    it2 = this->mons.find(it->first);
    if (it2 == this->mons.end())
      this->mons[it->first] = it->second;
    else
      this->mons[it->first] += it->second;
  }
  return (*this);
}

template<typename R, typename S>
PolynomialFp<R,S> & PolynomialFp<R,S>::operator-=(const PolynomialFp<R,S> & other)
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  
  for (it = other.mons.begin(); it != other.mons.end(); it++) {
    it2 = this->mons.find(it->first);
    if (it2 == this->mons.end())
      this->mons[it->first] = -it->second;
    else
      this->mons[it->first] -= it->second;
  }
  return (*this);
}

template<typename R, typename S>
PolynomialFp<R,S> & PolynomialFp<R,S>::operator*=(const PolynomialFp<R,S> & other)
{
  // Here we have no advantage doing it in place)
  (*this) = (*this)*other;
  return (*this);
}

template<typename R, typename S>
PolynomialFp<R,S> & PolynomialFp<R,S>::operator*=(const FpElement<R,S> & a)
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::iterator it;
  for (it = this->mons.begin(); it != this->mons.end(); it++) {
    it->second *= a;
  }
  return (*this);
}

template<typename R, typename S>
FpElement<R,S>
PolynomialFp<R,S>::evaluate(const std::vector<FpElement<R,S> > & vec) const
{
  FpElement<R,S> res(this->GF, 0);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->mons.begin(); i != this->mons.end(); i++) {
    FpElement<R,S> prod = i->second;
    std::multiset<size_t>::const_iterator j;
    for (j = i->first.begin(); j != i->first.end(); j++) {
#ifdef DEBUG
      assert((*j) < vec.size());
#endif
      prod *= vec[*j];
    }
    res += prod;
  }

  return res;
}

template<typename R, typename S>
PolynomialFp<R,S>
PolynomialFp<R,S>::evaluate(const std::vector<PolynomialFp<R,S> > & vec) const
{
  PolynomialFp<R,S> res(this->GF);

  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->mons.begin(); i != this->mons.end(); i++) {
    PolynomialFp<R,S> prod = i->second;
    std::multiset<size_t>::const_iterator j;
    for (j = i->first.begin(); j != i->first.end(); j++) {
#ifdef DEBUG
      assert((*j) < vec.size());
#endif
      prod *= vec[*j];
    }
    res += prod;
  }

  return res;
}

// booleans

template<typename R, typename S>
bool PolynomialFp<R,S>::is_zero() const
{
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  for (i = this->mons.begin(); i != this->mons.end(); i++)
    if (i->second != 0)
      return false;

  return true;
}

template<typename R, typename S>
bool PolynomialFp<R,S>::operator==(const PolynomialFp<R, S> & other) const
{
  return ((*this)-other).is_zero();
}

template<typename R, typename S>
bool PolynomialFp<R,S>::operator!=(const PolynomialFp<R, S> & other) const
{
  return !((*this)==other);
}

template<typename R, typename S>
bool PolynomialFp<R,S>::operator==(const FpElement<R, S> & a) const
{
  PolynomialFp<R,S> poly(a);
  return ((*this)==poly);
}

template<typename R, typename S>
bool PolynomialFp<R,S>::operator!=(const FpElement<R, S> & a) const
{
  return !((*this)==a);
}

template<typename R, typename S>
std::ostream& operator<<(std::ostream& os, const PolynomialFp<R,S>& poly)
{
  bool first = true;
  typename std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = poly.monomials().begin(); i != poly.monomials().end(); i++) {
    if (i->second == 0)
      continue;
    if (!first)
      os << "+";
    if ((i->second != 1) || (i->first).empty())
      os << i->second;
    std::multiset<size_t>::const_iterator j;
    bool inner_first = true;
    for (j = i->first.begin(); j != i->first.end(); j++) {
      if (!inner_first)
	os << "*";
      os << "x_" << (*j);
      inner_first = false;
    }
    first = false;
  }

  if (first)
    os << "0";

  return os;
}
