#include "birch_util.h"
#include <unordered_map>

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

// create the polynomial x^i
template<typename R>
UnivariatePoly<R> UnivariatePoly<R>::x(size_t i)
{
  UnivariatePoly<R> p;
  p.coeffs.resize(i+1);
  for (size_t j = 0; j < i; j++)
    p.coeffs[j] = Math<R>::zero();
  p.coeffs[i] = Math<R>::one();
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

template<typename R>
R UnivariatePoly<R>::content() const
{
  R c = Math<R>::zero();
  for (size_t i = 0; i < this->coeffs.size(); i++)
    c = Math<R>::gcd(c, this->coeffs[i]);

  return c;
}

// conversion, assignment operator
template<typename R>
template<typename T>
UnivariatePoly<R> &
UnivariatePoly<R>::operator=(const UnivariatePoly<T> & other)
{
  
  this->coeffs.resize(other.degree()+1);
  for (int i = 0; i <= other.degree(); i++) 
    this->coeffs[i] = birch_util::convert_Integer<T,R>(other.coefficient(i));
  
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
UnivariatePoly<R>
UnivariatePoly<R>::operator/(const UnivariatePoly<R> & other) const
{
  UnivariatePoly<R> q,r;
  div_rem((*this),other,q,r);

  return q;
}

template<typename R>
UnivariatePoly<R>
UnivariatePoly<R>::operator%(const UnivariatePoly<R> & other) const
{
  UnivariatePoly<R> q,r;
  div_rem((*this),other,q,r);

  return r;
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
UnivariatePoly<R> UnivariatePoly<R>::operator/(const R & a) const
{
  UnivariatePoly<R> prod;
  prod.coeffs.resize(this->coeffs.size());
  for (size_t i = 0; i < this->coeffs.size(); i++)
    prod.coeffs[i] = this->coeffs[i] / a;
  
  return prod;
}

template<typename R>
UnivariatePoly<R> UnivariatePoly<R>::operator%(const R & a) const
{
  UnivariatePoly<R> prod;
  prod.coeffs.resize(this->coeffs.size());
  for (size_t i = 0; i < this->coeffs.size(); i++)
    prod.coeffs[i] = this->coeffs[i] % a;
  
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
UnivariatePoly<R> &
UnivariatePoly<R>::operator/=(const UnivariatePoly<R> & other)
{
  (*this) = (*this)/other;
  return (*this);
}

template<typename R>
UnivariatePoly<R> &
UnivariatePoly<R>::operator%=(const UnivariatePoly<R> & other)
{
  (*this) = (*this)%other;
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
UnivariatePoly<R>& UnivariatePoly<R>::operator/=(const R & a)
{
  for (size_t i = 0; i < this->coeffs.size(); i++)
    this->coeffs[i] /= a;
  
  return (*this);
}

template<typename R>
UnivariatePoly<R>& UnivariatePoly<R>::operator%=(const R & a)
{
  for (size_t i = 0; i < this->coeffs.size(); i++)
    this->coeffs[i] %= a;
  
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
template<typename S, typename T>
UnivariatePolyFp<S, T>
UnivariatePoly<R>::mod(std::shared_ptr< const Fp<S, T> > GF) const
{
  UnivariatePolyFp<S, T> ret(GF);
  for (size_t i = 0; i < this->coeffs.size(); i++)
    // ret.coeffs.push_back(GF->mod(this->coeffs[i]));
    ret += (GF->mod(this->coeffs[i]))*UnivariatePolyFp<S, T>::x(GF,i);
  
  return ret;
}

template<typename R>
UnivariatePoly<R> UnivariatePoly<R>::derivative() const
{
  UnivariatePoly<R> f_prime;
  f_prime.coeffs.resize(this->degree());
  for (size_t i = 1; i < this->coeffs.size(); i++)
    f_prime.coeffs[i-1] = i * this->coeffs[i];
  
  return f_prime;
}

template<typename R>
void UnivariatePoly<R>::div_rem(const UnivariatePoly<R> & f,
				const UnivariatePoly<R> & g,
				UnivariatePoly<R> & q,
				UnivariatePoly<R> & r)
{
#ifdef DEBUG
  assert(g != 0);
#endif

  UnivariatePoly<R> t;
  
  q = 0;
  r = f;

  // we will use pseudo-remainder
  if (f.degree() >= g.degree()) {
    size_t d = f.degree() + 1 - g.degree();
    R g_d = g.lead();
    g_d = Math<R>::one();
    for (size_t i = 0; i < d; i++)
      g_d *= g.lead(); 
    //    R g_d = Math<R>::pow(g.lead(), d);
    r *= g_d;
  }

  while ((r != 0) && (r.degree() >= g.degree())) {
    R lc = r.lead() / g.lead();
    t = lc * x(r.degree()-g.degree());
    q += t;
    r -= t*g;
  }

  return;
}

template<typename R>
UnivariatePoly<R> UnivariatePoly<R>::gcd(const UnivariatePoly<R> & f,
					 const UnivariatePoly<R> & g)
{
  UnivariatePoly<R> q, r_minus, r, r_plus;
  r_minus = f / f.content();
  r = g / g.content();
  
  while (r != 0) {
    div_rem(r_minus, r, q, r_plus);
    r_minus = r;
    r = r_plus / r_plus.content();
  }
  
  return r_minus;
}

template<typename R>
UnivariatePoly<R> UnivariatePoly<R>::xgcd(const UnivariatePoly<R> & f,
					  const UnivariatePoly<R> & g,
					  UnivariatePoly<R> & s,
					  UnivariatePoly<R> & t)
{
  UnivariatePoly<R> q, r_minus, r, r_plus;
  UnivariatePoly<R> s_minus, s_plus, t_minus, t_plus;
  s = Math<R>::one();
  s_minus = Math<R>::zero();
  t = Math<R>::zero();
  t_minus = Math<R>::one();
  
  r_minus = f / f.content();
  r = g / g.content();
  
  while (r != Math<R>::zero()) {
    div_rem(r_minus, r, q, r_plus);
    
    R c = r_plus.content();
    W64 e = r_minus.degree()+1-r.degree();
    R a = Math<R>::pow(r.lead(), e);
    
    r_minus = r;
    r = r_plus / c;
    s_plus = (a*s_minus - q*s) / c;
    t_plus = (a*t_minus - q*t) / c;
    s_minus = s;
    t_minus = t;
    s = s_plus;
    t = t_plus;
  }

  // finalize
  s = s_minus;
  t = t_minus;
  
  return r_minus;
}

// We employ Yun's algorithm for this part
template<typename R>
std::vector< UnivariatePoly<R> >
UnivariatePoly<R>::squarefree_factor() const
{
  std::vector< UnivariatePoly<R> > fac;

  const UnivariatePoly<R> & f = *this;
  UnivariatePoly<R> f_prime = f.derivative();
  
  UnivariatePoly<R> a, b, c, d;
  b = f, c = f_prime, d = f_prime;
  
  while (b != Math<R>::one()) {
    a = gcd(b, d);
    fac.push_back(a);
    b /= a;
    c = d / a;
    d = c - b.derivative();
  }
  
  return fac;
}

template<typename R>
R UnivariatePoly<R>::landau_mignotte() const
{
  R d = this->degree() / 2;
  R B = Math<R>::binomial_coefficient(d-1, d/2-1);

  R norm = Math<R>::zero();
  for (size_t i = 0; i < this->coeffs.size(); i++)
    norm += this->coefficient(i)*this->coefficient(i);

  // we might need ceiling here
  norm = Math<R>::binomial_coefficient(d-1,d/2)*sqrt(norm);

  return B + norm;
}

// Here we asume f = prod(u) mod p^i
// and sum((prod(u)/u_i) v_i = 1 mod p^i
// Want to lift to the same but mod p^(i+1)

template<typename R>
template<typename S, typename T>
void UnivariatePoly<R>::hensel_step(std::vector<UnivariatePoly<R> > & u,
				    std::vector<UnivariatePoly<R> > & v,
				    std::shared_ptr< const Fp<S,T> > GF,
				    size_t i) const
{
  R p = GF->prime();
  R p_i = Math<R>::pow(p, i);
  UnivariatePoly<R> prod = Math<R>::one();
  for (size_t i = 0; i < u.size(); i++) {
    prod *= u[i];
  }
  UnivariatePoly<R> sum = Math<R>::zero();
  
#ifdef DEBUG
  assert(this->lead() % p != 0);
  assert((this->lead() - u[0].lead()) % p == 0);
  assert(u.size() == v.size());
  UnivariatePoly<R> sum2 = Math<R>::zero();
  for (size_t i = 0; i < u.size(); i++) {
    sum2 += (prod / u[i])*v[i];
    if (i > 0)
      assert(u[i].lead() == Math<R>::one());
    assert(v[i].degree() < u[i].degree());
  }
  assert( ((*this)-prod) % p_i == 0);
  assert( (sum2-Math<R>::one()) % p_i == 0);
#endif  

  // step 1 - lift the u_i
  
  u[0].lead() = this->lead();
  
  UnivariatePoly<R> t = ((*this) - prod) / p_i;
  UnivariatePoly<R> r;

  UnivariatePolyFp<S,T> t_p = t.mod(GF);
  UnivariatePolyFp<S,T> tv_bar(GF);
  UnivariatePolyFp<S,T> u_bar(GF);
  UnivariatePolyFp<S,T> q_bar(GF);
  UnivariatePolyFp<S,T> r_bar(GF);
  for (size_t i = 0; i < u.size(); i++) {
    u_bar = u[i].mod(GF);
    tv_bar = t_p*v[i].mod(GF);
    UnivariatePolyFp<S,T>::div_rem(tv_bar, u_bar, q_bar, r_bar);
    r = r_bar.lift();
    u[i] += p_i * r;
  }

  // step 2 - lift the v_i
  for (size_t i = 0; i < u.size(); i++) {
    sum += (prod / u[i])*v[i];
  }
  sum = - (sum - Math<R>::one()) / p_i ;
  UnivariatePolyFp<S,T> s_p = sum.mod(GF);

  for (size_t i = 0; i < u.size(); i++) {
    tv_bar = s_p*v[i].mod(GF);
    UnivariatePolyFp<S,T>::div_rem(tv_bar, u[i].mod(GF), q_bar, r_bar);
    r = r_bar.lift();
    v[i] += p_i * r;
  }
 
  return;
}

// here we lift all the way: f = prod(g) mod p
// and we lift to f = prod(g_lift) mod p^a

template<typename R>
template<typename S, typename T>
std::vector< UnivariatePoly<R> >
UnivariatePoly<R>::hensel_lift(const std::vector<UnivariatePolyFp<S, T> > & g,
			       size_t a) const
{
  R p = g[0].field()->prime();
  std::vector< UnivariatePoly<R> > u, v;
  std::vector< UnivariatePolyFp<S,T> > v_bar;
  UnivariatePolyFp<S,T> t(g[0].field());
  UnivariatePolyFp<S,T> prod = g[0];
  for (size_t i = 1; i < g.size(); i++) {
    UnivariatePolyFp<S,T>::xgcd(prod,g[i],v_bar[i],t);
    for (size_t j = 0; j < i; j++)
      v_bar[j] *= t;
  }
  v.resize(g.size());
  u.resize(g.size());
  for (size_t i = 0; i < g.size(); i++) {
    v[i] = v_bar[i].lift();
    u[i] = g[i].lift();
  }

  for (size_t i = 1; i < a; i++) {
    this->hensel_step(u, v, g[0].field(), i);
  }

  return u;
}

// create the polynomial x^i
template<typename R, typename S>
UnivariatePolyFp<R,S>
UnivariatePolyFp<R,S>::x(std::shared_ptr< const Fp<R,S> > GF,
			 size_t i)
{
  UnivariatePolyFp<R,S> p(GF);
  p.coeffs.resize(i+1);
  for (size_t j = 0; j < i; j++)
    p.coeffs[j] = Math< FpElement<R,S> >::zero();
  p.coeffs[i] = Math< FpElement<R,S> >::one();
  return p;
}

template<typename R, typename S>
UnivariatePoly<R> UnivariatePolyFp<R,S>::lift() const
{
  UnivariatePoly<R> ret(this->degree()+1);
  for (size_t i = 0; i < this->coeffs.size(); i++)
    // ret.coeffs[i] = this->coeffs[i].lift();
    ret += this->coeffs[i].lift()*UnivariatePoly<R>::x(i);

  return ret;
}

template<typename R, typename S>
UnivariatePolyFp<R,S>
UnivariatePolyFp<R,S>::pow_mod(size_t m, const UnivariatePolyFp<R,S> & f) const
{
  FpElement<R,S> one(GF_, 1); 
  UnivariatePolyFp<R,S> res(one);
  UnivariatePolyFp<R,S> q(GF_);
  UnivariatePolyFp<R,S> r(GF_);

  for (size_t i = 0; i < m; i++) {
    res = (*this)*res;
    div_rem(res, f, q, r);
    res = r;
  }
  
  return res;
}

template<typename R, typename S>
UnivariatePolyFp<R,S>::UnivariatePolyFp(const UnivariatePoly< FpElement<R,S> > & f)
{
  this->GF_ = f.const_coefficient().field();
  for (int i = 0; i <= f.degree(); i++)
    this->coeffs[i] = f.coefficient(i);
}

template<typename R, typename S>
UnivariatePolyFp<R,S> &
UnivariatePolyFp<R,S>::operator=(const UnivariatePoly< FpElement<R,S> > & f)
{
  this->GF_ = f.const_coefficient().field();
  for (int i = 0; i <= f.degree(); i++)
    this->coeffs[i] = f.coefficient(i);
  return (*this);
}

template<typename R, typename S>
std::vector< UnivariatePolyFp<R,S> >
UnivariatePolyFp<R,S>::cz_eq_deg_partial_factor(size_t r) const
{
  
  if (this->degree() == static_cast<int>(r)) {
    std::vector< UnivariatePolyFp<R,S> > ret(1, *this);
    return ret;
  }
  
  VectorFp<R,S,3> shifts(GF_);
  shifts[0] = -1;
  shifts[1] = 0;
  shifts[2] = -1;
  
  while (true) {
    UnivariatePolyFp<R,S> b(GF_);
    
    for (int i = 0; i < this->degree(); i++)
      b.coeffs.push_back(GF_->random());

    R p_r = Math<R>::pow(GF_->prime(),r);
    size_t m = (birch_util::convert_Integer<R, size_t>(p_r) - 1) / 2;

    UnivariatePolyFp<R,S> b_m = b.pow_mod(m, *this);
    UnivariatePolyFp<R,S> factor(GF_);
    UnivariatePolyFp<R,S> b_m_shifted(GF_);
    for (size_t i = 0; i < 3; i++) {
      b_m_shifted = b_m + shifts[i];
      factor = gcd(b_m_shifted, *this);
      if ((factor.degree() != 0) && (factor.degree() != this->degree()))
	return factor.cz_eq_deg_factor(r);
    }
  }
}

template<typename R, typename S>
std::vector< UnivariatePolyFp<R,S> >
UnivariatePolyFp<R,S>::cz_eq_deg_factor(size_t r) const
{
  std::vector< UnivariatePolyFp<R,S> > facs, partial;
  UnivariatePolyFp f = (*this);

  for (size_t i = 0; (r*i < this->coeffs.size()) && (f.degree() > 0); i++) {
    if (f.degree() == static_cast<int>(r)) {
      facs.push_back(f);
      return facs;
    }
    partial = f.cz_eq_deg_partial_factor(r);
    for ( UnivariatePolyFp<R,S> g : partial) {
      f /= g;
      facs.push_back(g);
    }
  }
  return facs;
}

template<typename R, typename S>
std::vector< UnivariatePolyFp<R,S> >
UnivariatePolyFp<R,S>::cz_distinct_deg_factor() const
{
  FpElement<R,S> one(GF_,Math<R>::one());
  size_t n = this->degree();
  std::vector< UnivariatePolyFp<R,S> > facs(n, one);
  
  
  size_t beta = n / 2;
  size_t l = floor(sqrt(beta));
  size_t m = (beta + l - 1) / l;
  R p = GF_->prime();
  R p_l = Math<R>::pow(p, l);

  std::vector< UnivariatePolyFp<R,S> > h, H, I;

  UnivariatePolyFp<R,S> x_p_i = UnivariatePolyFp<R,S>::x(GF_);
  for (size_t i = 0; i <= l+1; i++) {
    h.push_back(x_p_i);
    x_p_i = x_p_i.pow_mod(birch_util::convert_Integer<R, size_t>(p), *this);
  }

  x_p_i = x(GF_);
  for (size_t i = 0; i <= m+1; i++) {
    H.push_back(x_p_i);
    x_p_i = x_p_i.pow_mod(birch_util::convert_Integer<R, size_t>(p_l), *this);
  }

  UnivariatePolyFp<R,S> prod(GF_);
  UnivariatePolyFp<R,S> q(GF_);
  UnivariatePolyFp<R,S> r(GF_);
  UnivariatePolyFp<R,S> g(GF_);
  UnivariatePolyFp<R,S> mul(GF_);
  
  for (size_t i = 0; i <= m+1; i++) {
    prod = one;
    for (size_t j = 0; j < l; j++) {
      mul = prod*(H[i]-h[j]);
      div_rem(mul, *this, q, r);
      prod = r;
    }
    I.push_back(prod);
  }

  UnivariatePolyFp<R, S> f = *this;
  UnivariatePolyFp<R, S> diff(GF_);
  for (size_t i = 0; i <= m; i++) {
    g = gcd(f, I[i]);
    f /= g;
    for (size_t j = l; j > 0; j--) {
      diff = H[i] - h[j-1];
      facs[l*i-j] = UnivariatePolyFp<R,S>::gcd(g, diff);
      g /= facs[l*i-j];
    }
  }

  if (f.degree() >= 1)
    facs[f.degree()-1] = f;

  return facs;
}

template<typename R, typename S>
std::vector< UnivariatePolyFp<R,S> >
UnivariatePolyFp<R,S>::sqf_factor() const
{
  std::vector< UnivariatePolyFp<R,S> > fac;

  std::vector< UnivariatePolyFp<R,S> > dist = this->cz_distinct_deg_factor();

  for (size_t r = 0; r < dist.size(); r++) {
    std::vector< UnivariatePolyFp<R,S> > eq_deg =
      dist[r].cz_eq_deg_factor(r);

    fac.insert(fac.end(), eq_deg.begin(), eq_deg.end());
    
  }
  
  return fac;
}

template<typename R, typename S>
UnivariatePolyFp<R,S>
UnivariatePolyFp<R,S>::gcd(const UnivariatePolyFp<R,S> & f,
			   const UnivariatePolyFp<R,S> & g)
{
  const UnivariatePoly< FpElement<R,S> > & f1 =
    static_cast<const UnivariatePoly< FpElement<R,S> > &>(f);
  const UnivariatePoly< FpElement<R,S> > & g1 =
    static_cast<const UnivariatePoly< FpElement<R,S> > &>(g);
  
  UnivariatePoly< FpElement<R,S> > d = gcd(f1,g1);

  UnivariatePolyFp<R,S> res(d.const_coefficient().field());
  for (int i = 0; i <= d.degree(); i++)
    res.coeffs.push_back(d.coefficient(i));
  
  return res;
}

template<typename R, typename S>
UnivariatePolyFp<R,S>
UnivariatePolyFp<R,S>::xgcd(const UnivariatePolyFp<R,S> & f,
			    const UnivariatePolyFp<R,S> & g,
			    UnivariatePolyFp<R,S> & s,
			    UnivariatePolyFp<R,S> & t)
{
  const UnivariatePoly< FpElement<R,S> > & f1 =
    static_cast<const UnivariatePoly< FpElement<R,S> > &>(f);
  const UnivariatePoly< FpElement<R,S> > & g1 =
    static_cast<const UnivariatePoly< FpElement<R,S> > &>(g);
  /*
  UnivariatePoly< FpElement<R,S> > & s1 =
    static_cast< UnivariatePoly< FpElement<R,S> > &>(s);
  UnivariatePoly< FpElement<R,S> > & t1 =
    static_cast< UnivariatePoly< FpElement<R,S> > &>(t);
  
  UnivariatePoly< FpElement<R,S> > d = xgcd(f1,g1,s1,t1);

  UnivariatePolyFp<R,S> res(d.field());
  for (int i = 0; i <= d.degree(); i++)
    res.coeffs.push_back(d.coefficient(i));
  
  return res;
  */
  return xgcd(f1,g1,s,t);
}

template<typename R, typename S>
void UnivariatePolyFp<R,S>::div_rem(const UnivariatePolyFp<R,S> & f,
				    const UnivariatePolyFp<R,S> & g,
				    UnivariatePolyFp<R,S> & q,
				    UnivariatePolyFp<R,S> & r)
{
  const UnivariatePoly< FpElement<R,S> > & f1 =
    static_cast<const UnivariatePoly< FpElement<R,S> > &>(f);
  const UnivariatePoly< FpElement<R,S> > & g1 =
    static_cast<const UnivariatePoly< FpElement<R,S> > &>(g);
  /*
  UnivariatePoly< FpElement<R,S> > & q1 =
    static_cast< UnivariatePoly< FpElement<R,S> > &>(q);
  UnivariatePoly< FpElement<R,S> > & r1 =
    static_cast< UnivariatePoly< FpElement<R,S> > &>(r);
  
  div_rem(f1,g1,q1,r1);
  */
  div_rem(f1, g1, q, r);
  return;
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
