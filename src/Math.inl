// For future implementations that might be here

template <typename R>
Rational<R> Math<R>::bernoulli_number(const size_t & n)
{
  std::vector< Rational<R> > b = bernoulli_up_to(n);
  return b[n];
}

/* old code
template <typename R>
R Math<R>::binomial_coefficient(const R & n, const R & k)
{
  Rational<R> prod = Math<R>::one();
  for (R i = 0; i < k; i++)
    prod *= (n-i)/(k-i);
  return prod.floor();
}
*/

template<typename R>
R Math<R>::binomial_coefficient(const R & n, const R & k)
{
  R res = Math<R>::one();
  if (k > n - k)
    return binomial_coefficient(n, n-k);
  for (R i = 0; i < k; i++) {
    res *= (n-i);
    res /= (i+1);
  }
  return res;
}

template <typename R>
std::vector< Rational<R> > Math<R>::bernoulli_poly(const size_t & n)
{
  std::vector< Rational<R> > b = bernoulli_up_to(n);
  std::reverse(b.begin(), b.end());
  for (size_t k = 0; k <= n; k++)
    b[k] *= binomial_coefficient(n,k);
  return b;
}

// B_{n. chi} where chi is the quadratic character corresponding to
// quadratic field with discrminant d (Is it? verify we are working
// with the correct discriminant (d or 4d maybe?))

template <typename R>
int Math<R>::kronecker_symbol(const R & a, const R & n)
{
  // extremal cases
  if (n == 0) return (abs(a) == 1) ? 1 : 0;
  if (n == -1) return (a < 0) ? -1 : 1;
  if (n == 1) return 1;
  if (n == 2) {
    if (a % 2 == 0) return 0;
    R val = (a % 8)/2;
    if ((val == 0) || (val == 3))
      return 1;
    return -1;
  }
  if (a == -1) {
    R n_prime = n;
    while (n_prime % 2 == 0) n_prime /= 2;
    return ((n_prime / 2) % 2 == 0) ? 1 : -1;
  }
  if ((a == 2) && (n % 2 == 1)) {
    return ((n^2 / 8) % 2 == 0) ? 1 : -1;
  }
  // multiplicativity
  if (n < 0) return kronecker_symbol(a,-1)*kronecker_symbol(a,-n);
  if (a < 0) return kronecker_symbol(-1,n)*kronecker_symbol(-a,n);

  // now may assume n >= 3, a >= 0
 
  // quadratic reciprocity
  if (a < n) {
    R n_star;
    R n_prime = n;
    while (n_prime % 2 == 0) n_prime /= 2;
    n_star = ((n_prime / 2) % 2 == 0) ? n : -n;
    return kronecker_symbol(n_star, a);
  }

  // now we may also assume a ge n

  // if n = 2 mod 4, we can't reduce, use multiplicativity again
  if (n % 4 == 2) return kronecker_symbol(a, n/2)*kronecker_symbol(a,2);
  // now we can reduce
  return kronecker_symbol(a % n, n);
}

// B_{n. chi} where chi is the quadratic character corresponding to
// quadratic field with discrminant d (Is it? verify we are working
// with the correct discriminant (d or 4d maybe?))
template <typename R>
Rational<R> Math<R>::bernoulli_number(const size_t & n, const R & d)
{
  std::vector< Rational<R> > b = bernoulli_poly(n);
  R d_pow = 1;
  for (size_t k = 0; k < n; k++) d_pow *= d;
  Rational<R> b_chi = Math<R>::zero();
  for (R a = 0; a < d; a++)
    {
      int chi_a = kronecker_symbol(a, n);      
      R a_pow = 1;
      Rational<R> s = Math<R>::zero();
      for (size_t k = 0; k <= n; k++)
	{
	  d_pow /= d;
	  s += b[k]*a_pow*d_pow;
	  a_pow *= a;
	}
      s *= chi_a;
      b_chi += s;
    }
  return b_chi;
}

template<typename R>
R Math<R>::gcd(const R & a, const R & b)
{
  if (b < Math<R>::zero()) return gcd(a,-b);
  if (a < Math<R>::zero()) return gcd(-a,b);
  if (a < b) return gcd(b,a);
  if (b == Math<R>::zero()) return a;
  return gcd(b, a % b);
}

template<typename R>
R Math<R>::xgcd(const R & a, const R & b, R & x, R & y)
{
  R d;
  if (b < 0) {
    d = xgcd(a, -b, x, y);
    y = -y;
    return d;
  }
  if (a < 0) {
    d = xgcd(-a, b, x, y);
    x = -x;
    return d;
  }
  if (a < b) return xgcd(b, a, y, x);

  if (b == 0) {
    x = 1;
    y = 0;
    return a;
  }
  d = xgcd(b, a % b, y, x);
  y -= (a/b) * x;
  return d;
}

template<typename R>
R Math<R>::lcm(const R & a, const R & b)
{
  return a*b / gcd(a, b);
}

// Right now just a very naive implementation
template<typename R>
Rational<R> Math<R>::pow(const R & a, const Z64 & n)
{
  if (n < 0) return Math<R>::one()/pow(a,-n);
  Rational<R> ret = Math<R>::one();
  for (Z64 i = 0; i < n; i++) ret *= a;
  return ret;
}

template<typename R>
R Math<R>::pow(const R & a, const W64 & n)
{
  R ret = Math<R>::one();
  for (W64 i = 0; i < n; i++) ret *= a;
  return ret;
}
