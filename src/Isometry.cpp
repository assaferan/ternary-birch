#include "birch.h"
#include "birch_util.h"
#include "Isometry.h"

template<>
void Z_Isometry<3>::set_identity(void)
{
    mpz_set_ui(a11.get_mpz_t(), 1);
    mpz_set_ui(a12.get_mpz_t(), 0);
    mpz_set_ui(a13.get_mpz_t(), 0);
    mpz_set_ui(a21.get_mpz_t(), 0);
    mpz_set_ui(a22.get_mpz_t(), 1);
    mpz_set_ui(a23.get_mpz_t(), 0);
    mpz_set_ui(a31.get_mpz_t(), 0);
    mpz_set_ui(a32.get_mpz_t(), 0);
    mpz_set_ui(a33.get_mpz_t(), 1);
}

template<>
Z_Isometry<3> Z_Isometry<3>::inverse(const Z& p) const
{
  Z_Isometry<3> temp;

    mpz_mul(temp.a11.get_mpz_t(), a22.get_mpz_t(), a33.get_mpz_t());
    mpz_mul(temp.a12.get_mpz_t(), a13.get_mpz_t(), a32.get_mpz_t());
    mpz_mul(temp.a13.get_mpz_t(), a12.get_mpz_t(), a23.get_mpz_t());
    mpz_mul(temp.a21.get_mpz_t(), a23.get_mpz_t(), a31.get_mpz_t());
    mpz_mul(temp.a22.get_mpz_t(), a11.get_mpz_t(), a33.get_mpz_t());
    mpz_mul(temp.a23.get_mpz_t(), a13.get_mpz_t(), a21.get_mpz_t());
    mpz_mul(temp.a31.get_mpz_t(), a21.get_mpz_t(), a32.get_mpz_t());
    mpz_mul(temp.a32.get_mpz_t(), a12.get_mpz_t(), a31.get_mpz_t());
    mpz_mul(temp.a33.get_mpz_t(), a11.get_mpz_t(), a22.get_mpz_t());

    mpz_submul(temp.a11.get_mpz_t(), a23.get_mpz_t(), a32.get_mpz_t());
    mpz_submul(temp.a12.get_mpz_t(), a12.get_mpz_t(), a33.get_mpz_t());
    mpz_submul(temp.a13.get_mpz_t(), a13.get_mpz_t(), a22.get_mpz_t());
    mpz_submul(temp.a21.get_mpz_t(), a21.get_mpz_t(), a33.get_mpz_t());
    mpz_submul(temp.a22.get_mpz_t(), a13.get_mpz_t(), a31.get_mpz_t());
    mpz_submul(temp.a23.get_mpz_t(), a11.get_mpz_t(), a23.get_mpz_t());
    mpz_submul(temp.a31.get_mpz_t(), a22.get_mpz_t(), a31.get_mpz_t());
    mpz_submul(temp.a32.get_mpz_t(), a11.get_mpz_t(), a32.get_mpz_t());
    mpz_submul(temp.a33.get_mpz_t(), a12.get_mpz_t(), a21.get_mpz_t());

    mpz_divexact(temp.a11.get_mpz_t(), temp.a11.get_mpz_t(), p.get_mpz_t());
    mpz_divexact(temp.a12.get_mpz_t(), temp.a12.get_mpz_t(), p.get_mpz_t());
    mpz_divexact(temp.a13.get_mpz_t(), temp.a13.get_mpz_t(), p.get_mpz_t());
    mpz_divexact(temp.a21.get_mpz_t(), temp.a21.get_mpz_t(), p.get_mpz_t());
    mpz_divexact(temp.a22.get_mpz_t(), temp.a22.get_mpz_t(), p.get_mpz_t());
    mpz_divexact(temp.a23.get_mpz_t(), temp.a23.get_mpz_t(), p.get_mpz_t());
    mpz_divexact(temp.a31.get_mpz_t(), temp.a31.get_mpz_t(), p.get_mpz_t());
    mpz_divexact(temp.a32.get_mpz_t(), temp.a32.get_mpz_t(), p.get_mpz_t());
    mpz_divexact(temp.a33.get_mpz_t(), temp.a33.get_mpz_t(), p.get_mpz_t());

    return temp;
}

template<>
Z_Isometry<3> Z_Isometry<3>::operator*(const Z_Isometry<3>& s) const
{
   Z_Isometry<3> temp;

    mpz_mul(temp.a11.get_mpz_t(), a11.get_mpz_t(), s.a11.get_mpz_t());
    mpz_addmul(temp.a11.get_mpz_t(), a12.get_mpz_t(), s.a21.get_mpz_t());
    mpz_addmul(temp.a11.get_mpz_t(), a13.get_mpz_t(), s.a31.get_mpz_t());

    mpz_mul(temp.a12.get_mpz_t(), a11.get_mpz_t(), s.a12.get_mpz_t());
    mpz_addmul(temp.a12.get_mpz_t(), a12.get_mpz_t(), s.a22.get_mpz_t());
    mpz_addmul(temp.a12.get_mpz_t(), a13.get_mpz_t(), s.a32.get_mpz_t());

    mpz_mul(temp.a13.get_mpz_t(), a11.get_mpz_t(), s.a13.get_mpz_t());
    mpz_addmul(temp.a13.get_mpz_t(), a12.get_mpz_t(), s.a23.get_mpz_t());
    mpz_addmul(temp.a13.get_mpz_t(), a13.get_mpz_t(), s.a33.get_mpz_t());

    mpz_mul(temp.a21.get_mpz_t(), a21.get_mpz_t(), s.a11.get_mpz_t());
    mpz_addmul(temp.a21.get_mpz_t(), a22.get_mpz_t(), s.a21.get_mpz_t());
    mpz_addmul(temp.a21.get_mpz_t(), a23.get_mpz_t(), s.a31.get_mpz_t());

    mpz_mul(temp.a22.get_mpz_t(), a21.get_mpz_t(), s.a12.get_mpz_t());
    mpz_addmul(temp.a22.get_mpz_t(), a22.get_mpz_t(), s.a22.get_mpz_t());
    mpz_addmul(temp.a22.get_mpz_t(), a23.get_mpz_t(), s.a32.get_mpz_t());

    mpz_mul(temp.a23.get_mpz_t(), a21.get_mpz_t(), s.a13.get_mpz_t());
    mpz_addmul(temp.a23.get_mpz_t(), a22.get_mpz_t(), s.a23.get_mpz_t());
    mpz_addmul(temp.a23.get_mpz_t(), a23.get_mpz_t(), s.a33.get_mpz_t());

    mpz_mul(temp.a31.get_mpz_t(), a31.get_mpz_t(), s.a11.get_mpz_t());
    mpz_addmul(temp.a31.get_mpz_t(), a32.get_mpz_t(), s.a21.get_mpz_t());
    mpz_addmul(temp.a31.get_mpz_t(), a33.get_mpz_t(), s.a31.get_mpz_t());

    mpz_mul(temp.a32.get_mpz_t(), a31.get_mpz_t(), s.a12.get_mpz_t());
    mpz_addmul(temp.a32.get_mpz_t(), a32.get_mpz_t(), s.a22.get_mpz_t());
    mpz_addmul(temp.a32.get_mpz_t(), a33.get_mpz_t(), s.a32.get_mpz_t());

    mpz_mul(temp.a33.get_mpz_t(), a31.get_mpz_t(), s.a13.get_mpz_t());
    mpz_addmul(temp.a33.get_mpz_t(), a32.get_mpz_t(), s.a23.get_mpz_t());
    mpz_addmul(temp.a33.get_mpz_t(), a33.get_mpz_t(), s.a33.get_mpz_t());

    return temp;
}

template<>
void Z_Isometry<3>::A101011001()
{
    mpz_add(a13.get_mpz_t(), a13.get_mpz_t(), a12.get_mpz_t());
    mpz_add(a23.get_mpz_t(), a23.get_mpz_t(), a22.get_mpz_t());
    mpz_add(a33.get_mpz_t(), a33.get_mpz_t(), a32.get_mpz_t());
    mpz_add(a13.get_mpz_t(), a13.get_mpz_t(), a11.get_mpz_t());
    mpz_add(a23.get_mpz_t(), a23.get_mpz_t(), a21.get_mpz_t());
    mpz_add(a33.get_mpz_t(), a33.get_mpz_t(), a31.get_mpz_t());
}

template<>
void Z_Isometry<3>::A1t0010001(const Z& t)
{
    mpz_addmul(a12.get_mpz_t(), t.get_mpz_t(), a11.get_mpz_t());
    mpz_addmul(a22.get_mpz_t(), t.get_mpz_t(), a21.get_mpz_t());
    mpz_addmul(a32.get_mpz_t(), t.get_mpz_t(), a31.get_mpz_t());
}

template<>
void Z_Isometry<3>::A10001t001(const Z& t)
{
    mpz_addmul(a13.get_mpz_t(), t.get_mpz_t(), a12.get_mpz_t());
    mpz_addmul(a23.get_mpz_t(), t.get_mpz_t(), a22.get_mpz_t());
    mpz_addmul(a33.get_mpz_t(), t.get_mpz_t(), a32.get_mpz_t());
}

template<>
void Z_Isometry<3>::A10t010001(const Z& t)
{
    mpz_addmul(a13.get_mpz_t(), t.get_mpz_t(), a11.get_mpz_t());
    mpz_addmul(a23.get_mpz_t(), t.get_mpz_t(), a21.get_mpz_t());
    mpz_addmul(a33.get_mpz_t(), t.get_mpz_t(), a31.get_mpz_t());
}

template<>
void Z_Isometry<3>::A0n0n0000n()
{
    mpz_swap(a11.get_mpz_t(), a12.get_mpz_t());
    mpz_swap(a21.get_mpz_t(), a22.get_mpz_t());
    mpz_swap(a31.get_mpz_t(), a32.get_mpz_t());
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
    mpz_neg(a13.get_mpz_t(), a13.get_mpz_t());
    mpz_neg(a23.get_mpz_t(), a23.get_mpz_t());
    mpz_neg(a33.get_mpz_t(), a33.get_mpz_t());
}

template<>
void Z_Isometry<3>::An0000n0n0()
{
    mpz_swap(a12.get_mpz_t(), a13.get_mpz_t());
    mpz_swap(a22.get_mpz_t(), a23.get_mpz_t());
    mpz_swap(a32.get_mpz_t(), a33.get_mpz_t());
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
    mpz_neg(a13.get_mpz_t(), a13.get_mpz_t());
    mpz_neg(a23.get_mpz_t(), a23.get_mpz_t());
    mpz_neg(a33.get_mpz_t(), a33.get_mpz_t());
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
}

template<>
void Z_Isometry<3>::An00010001()
{
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
}

template<>
void Z_Isometry<3>::A1000n0001()
{
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
}

template<>
void Z_Isometry<3>::A10001000n()
{
    mpz_neg(a13.get_mpz_t(), a13.get_mpz_t());
    mpz_neg(a23.get_mpz_t(), a23.get_mpz_t());
    mpz_neg(a33.get_mpz_t(), a33.get_mpz_t());
}

template<>
void Z_Isometry<3>::An010n1001()
{
    mpz_add(a13.get_mpz_t(), a13.get_mpz_t(), a11.get_mpz_t());
    mpz_add(a13.get_mpz_t(), a13.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_add(a23.get_mpz_t(), a23.get_mpz_t(), a21.get_mpz_t());
    mpz_add(a23.get_mpz_t(), a23.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_add(a33.get_mpz_t(), a33.get_mpz_t(), a31.get_mpz_t());
    mpz_add(a33.get_mpz_t(), a33.get_mpz_t(), a32.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
}

template<>
void Z_Isometry<3>::Ann00n0001()
{
    mpz_add(a12.get_mpz_t(), a12.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_add(a22.get_mpz_t(), a22.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_add(a32.get_mpz_t(), a32.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
}

template<>
void Z_Isometry<3>::An0n01000n()
{
    mpz_add(a13.get_mpz_t(), a13.get_mpz_t(), a11.get_mpz_t());
    mpz_add(a23.get_mpz_t(), a23.get_mpz_t(), a21.get_mpz_t());
    mpz_add(a33.get_mpz_t(), a33.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a13.get_mpz_t(), a13.get_mpz_t());
    mpz_neg(a23.get_mpz_t(), a23.get_mpz_t());
    mpz_neg(a33.get_mpz_t(), a33.get_mpz_t());
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
}

template<>
void Z_Isometry<3>::A1000nn00n()
{
    mpz_add(a13.get_mpz_t(), a13.get_mpz_t(), a12.get_mpz_t());
    mpz_add(a23.get_mpz_t(), a23.get_mpz_t(), a22.get_mpz_t());
    mpz_add(a33.get_mpz_t(), a33.get_mpz_t(), a32.get_mpz_t());
    mpz_neg(a13.get_mpz_t(), a13.get_mpz_t());
    mpz_neg(a23.get_mpz_t(), a23.get_mpz_t());
    mpz_neg(a33.get_mpz_t(), a33.get_mpz_t());
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
}

template<>
void Z_Isometry<3>::Ann001000n()
{
    mpz_sub(a12.get_mpz_t(), a12.get_mpz_t(), a11.get_mpz_t());
    mpz_sub(a22.get_mpz_t(), a22.get_mpz_t(), a21.get_mpz_t());
    mpz_sub(a32.get_mpz_t(), a32.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a13.get_mpz_t(), a13.get_mpz_t());
    mpz_neg(a23.get_mpz_t(), a23.get_mpz_t());
    mpz_neg(a33.get_mpz_t(), a33.get_mpz_t());
}

template<>
void Z_Isometry<3>::An0n0n0001()
{
    mpz_sub(a13.get_mpz_t(), a13.get_mpz_t(), a11.get_mpz_t());
    mpz_sub(a23.get_mpz_t(), a23.get_mpz_t(), a21.get_mpz_t());
    mpz_sub(a33.get_mpz_t(), a33.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
}

template<>
void Z_Isometry<3>::An000nn001()
{
    mpz_sub(a13.get_mpz_t(), a13.get_mpz_t(), a12.get_mpz_t());
    mpz_sub(a23.get_mpz_t(), a23.get_mpz_t(), a22.get_mpz_t());
    mpz_sub(a33.get_mpz_t(), a33.get_mpz_t(), a32.get_mpz_t());
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
    mpz_neg(a11.get_mpz_t(), a11.get_mpz_t());
    mpz_neg(a21.get_mpz_t(), a21.get_mpz_t());
    mpz_neg(a31.get_mpz_t(), a31.get_mpz_t());
}

template<>
void Z_Isometry<3>::A1000010n0()
{
    mpz_swap(a13.get_mpz_t(), a12.get_mpz_t());
    mpz_swap(a23.get_mpz_t(), a22.get_mpz_t());
    mpz_swap(a33.get_mpz_t(), a32.get_mpz_t());
    mpz_neg(a12.get_mpz_t(), a12.get_mpz_t());
    mpz_neg(a22.get_mpz_t(), a22.get_mpz_t());
    mpz_neg(a32.get_mpz_t(), a32.get_mpz_t());
}

template<>
void Z_Isometry<3>::A1000100t1(const Z& t)
{
    mpz_addmul(a12.get_mpz_t(), a13.get_mpz_t(), t.get_mpz_t());
    mpz_addmul(a22.get_mpz_t(), a23.get_mpz_t(), t.get_mpz_t());
    mpz_addmul(a32.get_mpz_t(), a33.get_mpz_t(), t.get_mpz_t());
}

template<>
void Z_Isometry<3>::A100010t01(const Z& t)
{
    mpz_addmul(a11.get_mpz_t(), a13.get_mpz_t(), t.get_mpz_t());
    mpz_addmul(a21.get_mpz_t(), a23.get_mpz_t(), t.get_mpz_t());
    mpz_addmul(a31.get_mpz_t(), a33.get_mpz_t(), t.get_mpz_t());
}

template<>
void Z_Isometry<3>::A1000p000p2(const Z& p, const Z& pp)
{
    mpz_mul(a12.get_mpz_t(), a12.get_mpz_t(), p.get_mpz_t());
    mpz_mul(a22.get_mpz_t(), a22.get_mpz_t(), p.get_mpz_t());
    mpz_mul(a32.get_mpz_t(), a32.get_mpz_t(), p.get_mpz_t());
    mpz_mul(a13.get_mpz_t(), a13.get_mpz_t(), pp.get_mpz_t());
    mpz_mul(a23.get_mpz_t(), a23.get_mpz_t(), pp.get_mpz_t());
    mpz_mul(a33.get_mpz_t(), a33.get_mpz_t(), pp.get_mpz_t());
}
