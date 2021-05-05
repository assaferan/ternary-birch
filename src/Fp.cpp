#include "birch.h"
#include "Fp.h"

template<>
template<>
FpElement<W16, W32> W16_Fp::mod(const Z& a) const
{
    Z r;
    return FpElement(this, mpz_mod_ui(r.get_mpz_t(), a.get_mpz_t(), p));
}

template<>
template<>
FpElement<W32, W64> W32_Fp::mod(const Z& a) const
{
    Z r;
    return FpElement(this,mpz_mod_ui(r.get_mpz_t(), a.get_mpz_t(), p));
}

template<>
template<>
FpElement<W64, W128> W64_Fp::mod(const Z& a) const
{
    Z r;
    return FpElement(this,mpz_mod_ui(r.get_mpz_t(), a.get_mpz_t(), p));
}
