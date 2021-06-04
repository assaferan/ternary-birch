#include "Polynomial.h"

template<>
class UnivariatePolyFp<W16,W32>;

template<>
W64 UnivariatePoly<Z>::hash_value(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < this->coeffs.size(); i++)
    fnv = (fnv ^ mpz_get_si(this->coefficient(i).get_mpz_t())) * FNV_PRIME;
  return fnv;
}
