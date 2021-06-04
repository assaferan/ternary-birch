#include "Polynomial.h"

template class UnivariatePoly<Z>;
template class UnivariatePolyFp<W16,W32>;

template<>
W64 UnivariatePoly<Z>::hash_value(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < this->coeffs.size(); i++)
    fnv = (fnv ^ mpz_get_si(this->coefficient(i).get_mpz_t())) * FNV_PRIME;
  return fnv;
}

template<>
std::unordered_map< UnivariatePoly<Z>, size_t >
UnivariatePoly<Z>::factor() const
{
  std::unordered_map< UnivariatePoly<Z>, size_t > fac;

  std::vector< UnivariatePoly<Z> > sqf = this->squarefree_factor();
  std::random_device rd;
  W64 seed = rd();
  
  for (size_t i = 0; i < sqf.size(); i++) {
    UnivariatePoly<Z> f = sqf[i];
    if (f == Math<Z>::one()) continue;
    UnivariatePoly<Z> d = gcd(f, f.derivative());
    if (d == -1)
      d = 1;
    d -= Math<Z>::one();
    Z c = d.content();
    // for now we take an odd prime, to not have a special case
    // but in general, it might be bsest to work with 2
    Z p = Math<Z>::odd_prime_factor(c);
    W16 p_16 = birch_util::convert_Integer<Z, W16>(p);
    std::shared_ptr< const W16_Fp > GF = std::make_shared< W16_Fp >(p_16,seed);
    UnivariatePolyFp<W16,W32> f_p = f.mod(GF);
    std::vector< UnivariatePolyFp<W16, W32> > fac_p = f_p.sqf_factor();
    Z L = f.landau_mignotte();
    size_t a = 1;
    Z p_a = p;
    while (p_a <= 2*L) {
      a++;
      p_a *= p;
    }
    std::vector< UnivariatePoly<Z> > fac_lift = f.hensel_lift(fac_p,a);
    for ( UnivariatePoly<Z> g : fac_lift) {
      fac[g] = i;
    }
  }
  
  return fac;
}

template<>
W64 UnivariatePoly< W16_FpElement >::hash_value(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < this->coeffs.size(); i++)
    fnv = (fnv ^ this->coefficient(i).lift()) * FNV_PRIME;
  return fnv;
}

template<>
std::unordered_map< UnivariatePoly< W16_FpElement >, size_t >
UnivariatePoly< W16_FpElement >::factor() const
{
  // This is just a stub
  std::unordered_map< UnivariatePoly< W16_FpElement >, size_t > ret;
  return ret;
}
