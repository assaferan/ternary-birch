#ifndef __BIRCH_UTIL_H_
#define __BIRCH_UTIL_H_
// Birch-specific namespace

#include <map>

#include "birch.h"

namespace birch_util
{
  template<typename From, typename To>
  To convert_Integer(const From& x);

  template<typename From, typename To>
  PrimeSymbol<To> convert_PrimeSymbol(const PrimeSymbol<From>& symbol)
  {
    PrimeSymbol<To> temp;
    temp.p = convert_Integer<From,To>(symbol.p);
    temp.power = symbol.power;
    temp.ramified = symbol.ramified;
    return temp;
  }

  template<typename From, typename To, size_t n>
  QuadForm<To,n> convert_QuadForm(const QuadForm<From,n>& q)
  {
    SquareMatrix<To, n> mat;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	mat(i,j) = convert_Integer<From,To>(q.bilinear_form()(i,j));
    QuadForm<To,n> qq(mat);
    if (q.is_reduced()) {
      qq.set_reduced();
      qq.set_num_aut(qq.num_automorphisms());
    }
    
    return qq;
  }

  template<typename From, typename To, size_t n>
  Isometry<To,n> convert_Isometry(const Isometry<From,n>& s)
  {
    SquareMatrix<To, n> mat;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	mat(i,j) = convert_Integer<From,To>(s.a(i,j));
    To scale = convert_Integer<From, To>(s.get_scale());
    Isometry<To,n> ss(mat, scale);
    return ss;
  }

  template<typename From, typename To, size_t n>
  GenusRep<To,n> convert_GenusRep(const GenusRep<From,n>& from)
  {
    GenusRep<To,n> to;
    to.q = birch_util::convert_QuadForm<From,To,n>(from.q);
    to.s = birch_util::convert_Isometry<From,To,n>(from.s);
    to.sinv = birch_util::convert_Isometry<From,To,n>(from.sinv);
    to.parent = from.parent;
    to.p = convert_Integer<From,To>(from.p);
    for (auto& pair : from.es)
      {
	to.es[convert_Integer<From,To>(pair.first)] = pair.second;
      }
    return to;
  }

  template<typename R>
  inline R my_pow(const std::map<R,int>& pairs)
  {
    R x = 1;
    for (const auto& pair : pairs)
      {
	for (int k=0; k<pair.second; k++)
	  {
	    x *= pair.first;
	  }
      }
    return x;
  }

  int popcnt(Z64 x);

  extern int char_vals[256];

  inline int char_val(W64 x)
  {
    if (x <= 0xff) return char_vals[x];
    int value = 1;
    while (x)
      {
	value *= char_vals[x&0xff];
	x >>= 8;
      }
    return value;
  }

  // This method will return the integer divisor using at most 4 compares
  // if the divisor is equal to 0-7. Otherwise, it will compute the integer
  // divisor using a div instruction in addition to the 4 compares.
  // Experiments show that this is considerably faster than just blindly
  // issuing a div instruction for arbitrary input by a considerable amount.
  //
  // TODO: Determine whether this is faster than division after casting both
  // dividends to long double.
  template<typename R>
  inline R dumb_div(R num, R den)
  {
    if (num < den) return 0;

    if (num < 5*den)
      {
	if (num < 3*den)
	  {
	    if (num < 2*den) return 1;
	    return 2;
	  }
	else
	  {
	    if (num < 4*den) return 3;
	    return 4;
	  }
      }
    else
      {
	if (num < 7*den)
	  {
	    if (num < 6*den) return 5;
	    return 6;
	  }
	else
	  {
	    if (num < 8*den) return 7;
	  }
      }

    return num/den;
  }
}

#endif // __BIRCH_UTIL_H_
