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
PolynomialFp<R,S>::PolynomialFp(std::shared_ptr<const Fp<R,S>> GF, size_t i)
{
  this->GF = GF;
  std::multiset<size_t> singleton;
  singleton.insert(i);
  FpElement<R,S> one(GF,1);
  this->mons[singleton] = one;
}

template<typename R, typename S>
template<size_t n>
PolynomialFp<R,S>::PolynomialFp(const SquareMatrixFp<R, S, n> & q)
{
  this->GF = q.field();
  
  this->const_term = 0;
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
const FpElement<R, S> &
PolynomialFp<R,S>::coefficient(const std::multiset<size_t> & mon) const
{
  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  
  it = this->mons.find(mon);
  if (it == mons.end()) {
    FpElement<R,S> zero(this->GF,0);
    return zero;
  }
  return it->second;
}

template<typename R, typename S>
const FpElement<R, S> & PolynomialFp<R,S>::const_coefficient() const {
  std::multiset<size_t> empty_set;
  return this->coefficient(empty_set);
}

// coefficient of x_i
template<typename R, typename S>
const FpElement<R, S> & PolynomialFp<R,S>::coefficient(size_t i) const {
  std::multiset<size_t> singleton;
  singleton.insert(i);
  return this->coefficient(singleton);
}

// coefficient of x_i x_j
template<typename R, typename S>
const FpElement<R, S> & PolynomialFp<R,S>::coefficient(size_t i, size_t j) const
{
  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
  std::multiset<size_t> mon;
  mon.insert(i);
  mon.insert(j);
  return this->coefficient(mon);
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
  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
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
  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  for (it = this->mons.begin(); it != this->mons.end(); it++) {
    sum.mons[it->first] = it->second;
  }
  for (it = other.mons.begin(); it != other.end(); it++) {
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
  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it, it2;
  for (it = this->mons.begin(); it != this->mons.end(); it++) {
    diff.mons[it->first] = it->second;
  }
  for (it = other.mons.begin(); it != other.end(); it++) {
    it2 = diff.mons.find(it->first);
    if (it2 == sum.mons.end())
      diff.mons[it->first] = zero;
    diff.mons[it->first] -= it->second;
  }
  
  return diff;
}

template<typename R, typename S>
PolynomialFp<R,S> PolynomialFp<R,S>::operator*(const FpElement<R,S> & a) const
{
  PolynomialFp<R,S> prod(this->GF);

  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator it;
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

  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i,j,loc;
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
FpElement<R,S>
PolynomialFp<R,S>::evaluate(const std::vector<FpElement<R,S> > & vec) const
{
  FpElement<R,S> res(this->GF, 0);

  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->mons.begin(); i != this->mons.end(); i++) {
    FpElement<R,S> prod = i->second;
    std::multiset<size_t>::const_iterator j;
    for (j = i->first.begin(); i->first.end(); i++) {
      for (size_t k = 0; k < j->second; k++) {
#ifdef DEBUG
	assert(j->first < vec.size());
#endif
	prod *= vec[j->first];
      }
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

  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->mons.begin(); i != this->mons.end(); i++) {
    PolynomialFp<R,S> prod = i->second;
    std::multiset<size_t>::const_iterator j;
    for (j = i->first.begin(); i->first.end(); i++) {
      for (size_t k = 0; k < j->second; k++) {
#ifdef DEBUG
	assert(j->first < vec.size());
#endif
	prod *= vec[j->first];
      }
    }
    res += prod;
  }

  return res;
}

template<typename R, typename S>
PolynomialFp<R,S> PolynomialFp<R,S>::quadratic_part() const {
  PolynomialFp<R,S> quad(this->GF);

  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->mons.begin(); i != this->mons.end(); i++) {
    if ((i->first).size() == 2)
      quad.mons[i->first] = i->second;
  }
  
  return quad;
}

template<typename R, typename S>
std::vector< FpElement<R,S> > linear_part(size_t rank) const
{
  std::vector< FpElement<R,S> > linear;
  for (size_t i = 0; i < rank; i++)
    linear.push_back(this->coefficient(i));
  return linear;
}

template<typename R, typename S>
std::ostream& operator<<(std::ostream& os, const PolynomialFp<R,S>& poly)
{
  bool first = true;
  std::map<std::multiset<size_t>, FpElement<R,S> >::const_iterator i;
  
  for (i = this->mons.begin(); i != this->mons.end(); i++) {
    if (!first)
      os << "+";
    os << i->second;
    std::multiset<size_t>::const_iterator j;
    bool inner_first = true;
    for (j = i->first.begin(); j != i->first.end(); j++) {
      if (!inner_first)
	os << "*";
      os << "x_" << j->first;
      if (j->second != 1)
	os << "^" << j->second;
      inner_first = false;
    }
    first = false;
  }

  if (first)
    os << "0";

  return os;
}
