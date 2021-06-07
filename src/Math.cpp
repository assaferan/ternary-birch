#include "birch_util.h"
#include "Fp.h"
#include "Math.h"

template class Math<Z>;
template class Math<Z64>;
template class Math<Z128>;

template class Math<W16>;

const std::vector<int> hilbert_lut_odd = { 1, 1, 1, 1,
                                           1, 1,-1,-1,
                                           1,-1, 1,-1,
                                           1,-1,-1, 1 };

const int tab32[32] = {0,  9,  1, 10, 13, 21,  2, 29,
		       11, 14, 16, 18, 22, 25,  3, 30,
		       8, 12, 20, 28, 15, 17, 24,  7,
		       19, 27, 23,  6, 26,  5,  4, 31};


const int tab64[64] = { 63,  0, 58,  1, 59, 47, 53,  2,
			60, 39, 48, 27, 54, 33, 42,  3,
			61, 51, 37, 40, 49, 18, 28, 20,
			55, 30, 34, 11, 43, 14, 22,  4,
			62, 57, 46, 52, 38, 26, 32, 41,
			50, 36, 17, 19, 29, 10, 13, 21,
			56, 45, 25, 31, 35, 16,  9, 12,
			44, 24, 15,  8, 23,  7,  6,  5};

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

template<>
int Z128_Math::hilbert_symbol(Z128 a, Z128 b, const Z128& p)
{
  return Z_Math::hilbert_symbol(birch_util::convert_Integer<Z128,Z>(a),
				birch_util::convert_Integer<Z128,Z>(b),
				birch_util::convert_Integer<Z128,Z>(p));
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

// again we use the fact that we only consider small numbers
template <typename R>
R Math<R>::odd_prime_factor(const R & a)
{
#ifdef DEBUG
  assert( a != 2 );
#endif
  size_t p = 3;
  while (a % p != 0) p++;

  return p; 
}

template <>
Z Math<Z>::next_prime(const Z & a)
{
  Z p = a;
  mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
  return p; 
}

template <>
Z64 Math<Z64>::next_prime(const Z64 & a)
{
  return birch_util::convert_Integer<Z, Z64>(Z_Math::next_prime(a));
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

template<>
Z Math<Z>::zero()
{
  return 0;
}

template<>
bool Math<Z>::is_zero(const Z & a)
{
  return (a == 0);
}

template<>
Z32 Math<Z32>::zero()
{
  return 0;
}

template<>
bool Math<Z32>::is_zero(const Z32 & a)
{
  return (a == 0);
}

template<>
W32 Math<W32>::zero()
{
  return 0;
}

template<>
bool Math<W32>::is_zero(const W32 & a)
{
  return (a == 0);
}
/*
template<>
NumberFieldElement<Z> Math< NumberFieldElement<Z> >::zero()
{
  return 0;
  }
*/

template<>
W16_FpElement Math< W16_FpElement >::zero()
{
  return 0;
}

template<>
bool Math<W16_FpElement>::is_zero(const W16_FpElement & a)
{
  return a.is_zero();
}

template<>
W32_FpElement Math< W32_FpElement >::zero()
{
  return 0;
}

template<>
bool Math<W32_FpElement>::is_zero(const W32_FpElement & a)
{
  return a.is_zero();
}


template<>
W64_FpElement Math< W64_FpElement >::zero()
{
  return 0;
}

template<>
bool Math<W64_FpElement>::is_zero(const W64_FpElement & a)
{
  return a.is_zero();
}

template<>
FpElement<W16, W16> Math< FpElement<W16, W16> >::zero()
{
  return 0;
}

template<>
bool Math<FpElement<W16,W16> >::is_zero(const FpElement<W16,W16> & a)
{
  return a.is_zero();
}


template<>
Z64 Math<Z64>::zero()
{
  return 0;
}

template<>
bool Math<Z64>::is_zero(const Z64 & a)
{
  return (a == 0);
}

template<>
Z128 Math<Z128>::zero()
{
  return 0;
}

template<>
bool Math<Z128>::is_zero(const Z128 & a)
{
  return (a == 0);
}

template<>
W16 Math<W16>::zero()
{
  return 0;
}

template<>
bool Math<W16>::is_zero(const W16 & a)
{
  return (a == 0);
}

template<>
W64 Math<W64>::zero()
{
  return 0;
}

template<>
bool Math<W64>::is_zero(const W64 & a)
{
  return (a == 0);
}

template<>
Rational<Z> Math< Rational<Z> >::zero()
{
  Z zero = 0;
  return zero;
}

template<>
bool Math<Rational<Z> >::is_zero(const Rational<Z> & a)
{
  return (a.num() == 0);
}

template<>
Rational<Z64> Math< Rational<Z64> >::zero()
{
  Z64 zero = 0;
  return zero;
}

template<>
bool Math<Rational<Z64> >::is_zero(const Rational<Z64> & a)
{
  return (a.num() == 0);
}

template<>
Rational<Z128> Math< Rational<Z128> >::zero()
{
  Z128 zero = 0;
  return zero;
}

template<>
bool Math<Rational<Z128> >::is_zero(const Rational<Z128> & a)
{
  return (a.num() == 0);
}

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

template<>
Z Math<Z>::one()
{
  return 1;
}

template<>
Z32 Math<Z32>::one()
{
  return 1;
}

template<>
Z64 Math<Z64>::one()
{
  return 1;
}

template<>
Z128 Math<Z128>::one()
{
  return 1;
}

template<>
W16 Math<W16>::one()
{
  return 1;
}

template<>
W32 Math<W32>::one()
{
  return 0;
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

template<>
Rational<Z128> Math< Rational<Z128> >::one()
{
  Z128 one = 1;
  return one;
}

template<>
int Math<uint64_t>::log2(const uint64_t & a)
{
  uint64_t value = a;
  value |= value >> 1;
  value |= value >> 2;
  value |= value >> 4;
  value |= value >> 8;
  value |= value >> 16;
  value |= value >> 32;
  return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}

template<>
int Math<uint32_t>::log2(const uint32_t & a)
{
  uint32_t value = a;
  value |= value >> 1;
  value |= value >> 2;
  value |= value >> 4;
  value |= value >> 8;
  value |= value >> 16;
  return tab32[(uint32_t)(value*0x07C4ACDD) >> 27];
}

template<>
int Math<long>::log2(const long & a)
{
#ifdef DEBUG
  assert(a > 0);
#endif
  unsigned long value = a;
  return Math<unsigned long>::log2(value);
}

template<>
int Math<Z>::log2(const Z & value)
{
  long exp;
  mpz_get_d_2exp(&exp, value.get_mpz_t());
  return exp;
}
