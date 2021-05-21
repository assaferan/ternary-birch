#include "Fp.h"
#include "Math.h"

template class Math<Z>;
template class Math<Z64>;

const std::vector<int> hilbert_lut_odd = { 1, 1, 1, 1,
                                           1, 1,-1,-1,
                                           1,-1, 1,-1,
                                           1,-1,-1, 1 };
// This table seems to be incorrect,
// Maybe the indexing is wrong.
// In any case, here is my attempt at a correct table
// In retrospect, there was one wrong value, which we seemed to hit :-)
/*
const std::vector<int> hilbert_lut_p2 = { 1, 1, 1, 1, 1, 1, 1, 1,
                                          1,-1, 1,-1,-1, 1,-1, 1,
                                          1, 1, 1, 1,-1,-1,-1,-1,
                                          1,-1, 1,-1, 1,-1, 1,-1,
                                          1,-1,-1, 1, 1,-1,-1, 1,
                                          1, 1,-1,-1,-1,-1,-1, 1,
                                          1,-1,-1, 1,-1, 1, 1,-1,
                                          1, 1,-1,-1, 1, 1,-1,-1 };
*/
const std::vector<int> hilbert_lut_p2 = { 1, 1, 1, 1, 1, 1, 1, 1,
                                          1,-1, 1,-1,-1, 1,-1, 1,
                                          1, 1, 1, 1,-1,-1,-1,-1,
                                          1,-1, 1,-1, 1,-1, 1,-1,
                                          1,-1,-1, 1, 1,-1,-1, 1,
                                          1, 1,-1,-1,-1,-1, 1, 1,
                                          1,-1,-1, 1,-1, 1, 1,-1,
                                          1, 1,-1,-1, 1, 1,-1,-1 };

template<>
int Z_Math::hilbert_symbol(Z a, Z b, const Z& p)
{
    int a_val = 0;
    while (a % p == 0)
    {
        ++a_val;
        a /= p;
    }

    int b_val = 0;
    while (b % p == 0)
    {
        ++b_val;
        b /= p;
    }

    if (p == 2)
    {
        int aa = (mpz_class(a%8).get_si() >> 1) & 0x3;
        int bb = (mpz_class(b%8).get_si() >> 1) & 0x3;
        int index = ((a_val&0x1)<<5) | (aa << 3) | ((b_val&0x1)<<2) | bb;
        return hilbert_lut_p2[index];
    }

    int a_notsqr = mpz_legendre(a.get_mpz_t(), p.get_mpz_t()) == -1;
    int b_notsqr = mpz_legendre(b.get_mpz_t(), p.get_mpz_t()) == -1;

    int index = ((a_val&0x1)<<3) | (a_notsqr<<2) | ((b_val&0x1)<<1) | b_notsqr;
    if (((index & 0xa) == 0xa) && ((p%4) == 0x3))
    {
        return -hilbert_lut_odd[index];
    }
    else
    {
        return hilbert_lut_odd[index];
    }
}

template<>
int Z64_Math::hilbert_symbol(Z64 a, Z64 b, const Z64& p)
{
    return Z_Math::hilbert_symbol(a, b, p);
}

// This is a helper function
// !! TODO - maybe have a utils file to put it there
// We do only trial divison, since our numbers are always small enough
// (num will be in the order of 1000 at most)
template <typename R>
std::vector< std::pair<R, size_t> > Math<R>::factorization(const R & num)
{
  std::vector< std::pair<R, size_t> > factors;
  R temp_num = abs(num);
  size_t exp;
  R a = 2;
  while (temp_num != 1)
    {
      if (temp_num % a == 0)
	{
	  exp = 1; 
	  temp_num /= a;
	  while (temp_num % a == 0)
	    {
	      exp++;
	      temp_num /= a;
	    }
	  factors.push_back(std::make_pair(a, exp));
	}
      a++;
    }
  return factors;
}


template <typename R>
size_t Math<R>::valuation(const R & a, const R& p)
{
  assert(a != 0);
  R t = a;
  
  size_t exp = 0;

  while (t % p == 0)
    {
      exp++;
      t /= p;
    }

  return exp;
}

template <typename R>
int Math<R>::valuation(const Rational<R> & a, const R& p)
{
  return valuation(a.num(), p) - valuation(a.denom(), p);
}

template <typename R>
bool Math<R>::is_local_square(const R & a, const R& p)
{
  if (a == 0) return true;
  size_t v = valuation(a,p);
  if (v % 2 == 1) return false;
  R a0 = a;
  for (size_t i = 0; i < v; i++) a0 /= p;
  bool ok = (kronecker_symbol(a0,p) != -1);
  if (p != 2) return ok;
  size_t w = valuation(a0-1,p);
  assert(w >= 1);
  size_t ee = 2;

  while ((w < ee) && (w % 2 == 0)) {
    a0 /= (1+ (1<< (w/2)))*(1+ (1<< (w/2)));
    w = valuation(a0-1,p);
  }
  return ((w > ee) || ((w == ee) && (a0 % 8 == 1)));
}

template <typename R>
bool Math<R>::is_local_square(const Rational<R> & a, const R& p)
{
  return (is_local_square(a.num(),p) == is_local_square(a.denom(),p));
}

template <typename R>
bool Math<R>::is_square(const R & num)
{
  // We could do that without factoring, but it's good enough for now
  std::vector< std::pair<R, size_t> > fac = factorization(num);
  if (num < 0) return false;
  for (std::pair<R, size_t> fa : fac)
    if (fa.second % 2 == 1) return false;
  return true;
}


// collect all of the along the way.
// We are using Akiyama and Tanigawa's algorithm
// It's not the fastest, but it is one of the simplest.
template <typename R>
std::vector< Rational<R> > Math<R>::bernoulli_up_to(const size_t & n)
{
  std::vector< Rational<R> > a(n+1);
  std::vector< Rational<R> > b(n+1);
  for (size_t i = 0; i <= n; i++) {
    Rational<R> r(1, i+1);
    a[i] = r;
  }
  b[0] = a[0];
  for (size_t i = 0; i < n; i++)
    {
      for (size_t j = 0; j < n - i; j++) {
	R mult = j + 1;
	a[j] = mult*(a[j] - a[j+1]);
      }
      b[i+1] = a[0];
    }
  
  return b;
  
}

template <typename R>
Rational<R> Math<R>::bernoulli_number(const size_t & n)
{
  std::vector< Rational<R> > b = bernoulli_up_to(n);
  return b[n];
}

template <typename R>
R Math<R>::binomial_coefficient(const R & n, const R & k)
{
  Rational<R> prod = Math<R>::one();
  for (R i = 0; i < k; i++)
    prod *= (n-i)/(k-i);
  return prod.floor();
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
  if (a < 0) return kronecker_symbol(-1,n)*kronecker_symbol(a,n);

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
  if (b < 0) return gcd(a,-b);
  if (a < 0) return gcd(-a,b);
  if (a < b) return gcd(b,a);
  if (b == 0) return a;
  return gcd(b, a % b);
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

// these seem to be covered by birch_util
/*
template<>
Z64 Math<Z>::get_int(const Z & a)
{
  return mpz_get_si(a.get_mpz_t());
}

template<>
Z64 Math<Z64>::get_int(const Z64 & a)
{
  return a;
}
*/

template<>
Z Math<Z>::zero()
{
  return 0;
}

/*
template<>
W16_FpElement Math< W16_FpElement >::zero()
{
  return 0;
}

template<>
W32_FpElement Math< W32_FpElement >::zero()
{
  return 0;
}

template<>
W64_FpElement Math< W64_FpElement >::zero()
{
  return 0;
}

template<>
FpElement<W16, W16> Math< FpElement<W16, W16> >::zero()
{
  return 0;
}
*/

template<>
Z64 Math<Z64>::zero()
{
  return 0;
}

template<>
W16 Math<W16>::zero()
{
  return 0;
}

template<>
W64 Math<W64>::zero()
{
  return 0;
}

template<>
Rational<Z> Math< Rational<Z> >::zero()
{
  Z zero = 0;
  return zero;
}

template<>
Rational<Z64> Math< Rational<Z64> >::zero()
{
  Z64 zero = 0;
  return zero;
}

/*
template<>
W16_FpElement Math< W16_FpElement >::one()
{
  return 1;
}

template<>
W32_FpElement Math< W32_FpElement >::one()
{
  return 1;
}

template<>
W64_FpElement Math< W64_FpElement >::one()
{
  return 1;
}

template<>
FpElement<W16, W16> Math< FpElement<W16, W16> >::one()
{
  return 1;
}
*/

template<>
Z Math<Z>::one()
{
  return 1;
}

template<>
Z64 Math<Z64>::one()
{
  return 1;
}

template<>
W16 Math<W16>::one()
{
  return 1;
}


template<>
W64 Math<W64>::one()
{
  return 1;
}

template<>
Rational<Z> Math< Rational<Z> >::one()
{
  Z one = 1;
  return one;
}

template<>
Rational<Z64> Math< Rational<Z64> >::one()
{
  Z64 one = 1;
  return one;
}

