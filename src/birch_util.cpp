#include "birch.h"
#include "birch_util.h"

namespace birch_util
{
  static int bitcounts[256] = {
			       0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2,
			       3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3,
			       3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3,
			       4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4,
			       3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5,
			       6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4,
			       4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
			       6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5,
			       3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3,
			       4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6,
			       6, 7, 6, 7, 7, 8
  };

  int popcnt(Z64 x)
  {
    if (x <= 0xff) return bitcounts[x];

    int count = 0;
    while (x)
      {
	count += bitcounts[x&0xff];
	x >>= 8;
      }
    return count;
  }

  template<>
  W16 convert_Integer<Z>(const Z& x)
  {
    return (W16)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W16 convert_Integer<Z64>(const Z64& x)
  {
    return (W16)x;
  }

  template<>
  W16 convert_Integer<Z128>(const Z128& x)
  {
    return (W16)x;
  }
  
  template<>
  W32 convert_Integer<Z>(const Z& x)
  {
    return (W32)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W64 convert_Integer<W16>(const W16& x)
  {
    return (W16)x;
  }
  
  template<>
  W32 convert_Integer<Z64>(const Z64& x)
  {
    return (W32)x;
  }

  template<>
  W64 convert_Integer<Z>(const Z& x)
  {
    return (W64)mpz_get_ui(x.get_mpz_t());
  }

  template<>
  W64 convert_Integer<Z64>(const Z64& x)
  {
    return (W64)x;
  }

  template<>
  Z convert_Integer<W16>(const W16& x)
  {
    return Z(x);
  }
  
  template<>
  Z convert_Integer<Z32>(const Z32& x)
  {
    return Z(x);
  }
  
  template<>
  Z convert_Integer<Z64>(const Z64& x)
  {
    return Z(x);
  }

  template<>
  Z convert_Integer<Z128>(const Z128& x)
  {
    Z64 top = x >> 64;
    Z64 bottom = x & UINT64_MAX;
    // !! TODO - check signs
    Z ret = top;
    ret <<= 64;
    ret += bottom;
    return ret;
  }

  template<>
  Z64 convert_Integer<Z>(const Z& x)
  {
    return mpz_get_si(x.get_mpz_t());
  }

  template<>
  Z64 convert_Integer<Z64>(const Z64& x)
  {
    return x;
  }

  template<>
  Z convert_Integer<Z>(const Z& x)
  {
    return x;
  }

  int char_vals[256] = {
			1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
			1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
			1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
			1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,
			1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,
			1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
			1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,
			1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,
			1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1,
			-1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,
			1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
			1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,
			1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1,
			-1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
			1, -1, -1,  1
  };
}
