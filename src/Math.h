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
  static int hilbert_symbol(R a, R b, const R& p);
  static int kronecker_symbol(const R & a, const R & n);
  static size_t valuation(const R& a, const R& p);
  static bool is_local_square(const R& a, const R& p);
  static std::vector< std::pair<R, size_t> > factorization(const R & num);
  static bool is_square(const R & num);
  static Rational<R> bernoulli_number(const size_t &);
  static R binomial_coefficient(const R & n, const R & k);
  static std::vector< Rational<R> > bernoulli_poly(const size_t & n);
  static Rational<R> bernoulli_number(const size_t & n, const R & d);
};

#include "Math.inl"

#endif // __MATH_H_
