#ifndef __MATH_H_
#define __MATH_H_

#include <vector>

// Is this one really needed?
#include "birch.h"
#include "Rational.h"

template<typename R>
class Math
{
public:
  static Rational<R> pow(const R & a, const R & n);
  static R gcd(const R & a, const R & b);
  static R lcm(const R & a, const R & b);
  static int hilbert_symbol(R a, R b, const R& p);
  static int kronecker_symbol(const R & a, const R & n);
  static size_t valuation(const R& a, const R& p);
  static int valuation(const Rational<R>& a, const R& p);
  static bool is_local_square(const R& a, const R& p);
  static bool is_local_square(const Rational<R>& a, const R& p);
  static std::vector< std::pair<R, size_t> > factorization(const R & num);
  static bool is_square(const R & num);
  static std::vector< Rational<R> > bernoulli_up_to(const size_t &);
  static Rational<R> bernoulli_number(const size_t &);
  static R binomial_coefficient(const R & n, const R & k);
  static std::vector< Rational<R> > bernoulli_poly(const size_t & n);
  static Rational<R> bernoulli_number(const size_t & n, const R & d);

  // static Z64 get_int(const R & a);

  static R zero();
  static R one();
};

#include "Math.inl"

#endif // __MATH_H_
