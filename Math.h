#ifndef __MATH_H_
#define __MATH_H_

#include <cmath>
#include <gmpxx.h>

template<typename R, typename F>
class Math
{
public:
    static R gcd(const R& a, const R& b);
    static R pow(const R& a, const R& e, const R& p);
    static R modinv(const R& a, const R& p);
    static void fix_vector(std::vector<R>& vec, const R& p);
    static void normalize_vector(std::vector<R>& vec, const R& p);
    static int64_t kronecker(const R& a, const R& b);
    static R square_root(R a, const R& p);
    static int64_t valuation(F x, R p);
    static std::vector<R> prime_divisors_naive(R n);
    static std::vector<R> prime_divisors(const R& n);
    static std::vector<std::vector<R>> squarefree_divisors(const R& n);
    static std::vector<R> prime_up_to(const R& upTo, const R& coprimeTo);
};

#endif // __MATH_H_
