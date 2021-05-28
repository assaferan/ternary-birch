#ifndef __RATIONAL_H_
#define __RATIONAL_H_

#include "Math.h"

template<typename R>
class Rational
{
public:
  
  // c-tors
  Rational(const R& num, const R& denom) : num_(num), denom_(denom)
  { reduce(); }

  Rational(const R & num) : num_(num), denom_(1) {}
  
  // default c-tor
  Rational() : num_(0), denom_(1) {}

  // copy c-tor
  Rational(const Rational<R> & other)
    : num_(other.num_), denom_(other.denom_) {}

  // access
  const R & num() const {return this->num_; }
  const R & denom() const {return this->denom_; }

  // arithmetic
  Rational<R> operator+() const {return Rational(num_, denom_); }
  Rational<R> operator-() const {return Rational(-num_, denom_); }
  Rational<R> operator+(const Rational<R> &) const;
  Rational<R> operator-(const Rational<R> &b) const {return (*this)+(-b); }
  Rational<R> operator*(const Rational<R> &) const;
  Rational<R> operator*(const R & b) const {
    Rational<R> b_rat(b);
    return (*this)*b_rat;
  }

  Rational<R> operator/(const Rational<R> &) const;
  Rational<R> operator/(const R & b) const {
    Rational<R> b_rat(b);
    return (*this)/b_rat;
  }
    
  // assignment
  Rational<R> & operator=(const Rational<R> & b)
  {
    if (this != &b) {
      num_ = b.num_;
      denom_ = b.denom_;
    }
    return (*this);
  }
  
  Rational<R> & operator+=(const Rational<R> &b)
  {return ((*this) = (*this) + b);}
  Rational<R> & operator-=(const Rational<R> &b)
  {return ((*this) = (*this) - b);}
  Rational<R> & operator*=(const Rational<R> &b)
  {return ((*this) = (*this) * b);}
  Rational<R> & operator*=(const R &b)
  {return ((*this) = (*this) * b);}
  Rational<R> & operator/=(const Rational<R> &b)
  {return ((*this) = (*this) / b);}
  Rational<R> & operator/=(const R &b)
  {return ((*this) = (*this) / b);}

  // comparison
  bool operator==(const Rational<R> &) const;
  bool operator!=(const Rational<R> &b) const {return !((*this)==b); }
  bool operator<(const Rational<R> &) const;
  bool operator>(const Rational<R> &b) const {return b < (*this); }
  bool operator<=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) < b); }
  bool operator>=(const Rational<R> &b) const
  {return ((*this) == b) || ((*this) > b); }

  // other
  R floor() const    
  {
    R num = num_;
    R denom = denom_;
    
    if (denom < 0) {
      denom = -denom;
      num = - num;
    }
      
    return ((num >= Math<R>::zero()) ? num : (num - denom + 1)) / denom;
  }

  R ceiling() const
  {
    R num = num_;
    R denom = denom_;
    
    if (denom < 0) {
      denom = -denom;
      num = - num;
    }
      
    return ((num >= Math<R>::zero()) ? (num + denom - 1) : num) / denom;

  }

  bool is_integral() const
  {R one = 1; return ((denom_ == one) || (denom_ == -one)); }

protected:
  R num_;
  R denom_;

private:
  void reduce(void);
};

// other
template <typename R>
Rational<R> operator*(R b, const Rational<R> & r) {
  return r*b;
}

template <typename R>
Rational<R> operator-(R b, const Rational<R> & r) {
  Rational<R> b_rat(b);
  return b_rat-r;
}

template <typename R>
Rational<R> operator+(R b, const Rational<R> & r) {
  Rational<R> b_rat(b);
  return b_rat+r;
}

template <typename R>
Rational<R> operator/(R b, const Rational<R> & r) {
  Rational<R> b_rat(b);
  return b_rat/r;
}

template<typename R>
static Rational<R> abs(const Rational<R> & r)
{ R zero = 0; return (r > zero) ? r : -r;}

template <typename R>
std::ostream& operator<<(std::ostream & os, const Rational<R> & r)
{
  R one = 1;
  if (r.denom() == one) return os << r.num();
  if (r.denom() == -one) return os << -r.num();
  os << r.num() << "/" << r.denom();
  return os;
}

#include "Rational.inl"

#endif // __RATIONAL_H_
