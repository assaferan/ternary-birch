#ifndef __SPINOR_H_
#define __SPINOR_H_

#include "birch.h"

template<typename R>
class Spinor
{
  template<typename T>
    friend class Spinor;

public:
    Spinor(const std::vector<R>& primes)
    {
        this->primes_ = primes;
        this->twist = (1LL << this->primes_.size()) - 1;
    }

  template<size_t n>
  Z64 norm(const QuadForm<R, n>& q, const Isometry<R, n>& s, const R& scalar) const
    {
      R tr = 0;
      for (size_t i = 0; i < n; i++)
	tr += s(i,i);
      // Stub
      // !! TODO - understand what should be happening here...
      // It seems it should feed compute_vals the spinor norm of s,
      // but there are subtleties
      // Also we would probably prefer to compute the character directly
      return this->compute_vals(tr + scalar);
    }

    const std::vector<R> primes(void) const
    {
        return this->primes_;
    }

private:
    std::vector<R> primes_;
    Z64 twist;

    Z64 compute_vals(R x) const
    {
        Z64 val = 0;
        Z64 mask = 1;
        for (const R& p : this->primes_)
        {
            while (x % p == 0)
            {
                x /= p;
                val ^= mask;
            }
            mask <<= 1;
        }
        return val;
    }
};

#endif // __SPINOR_H_
