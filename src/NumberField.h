#ifndef __NUMBER_FIELD_H_
#define __NUMBER_FIELD_H_

#include "Polynomial.h"
#include "Rational.h"

template<typename R>
class NumberFieldElement
{
public:
  
  NumberFieldElement(const UnivariatePoly< Rational<R> > & poly) : elt(poly) {}
  NumberFieldElement(const R & a) : elt(a) {} 
  NumberFieldElement(const Rational<R> & a) : elt(a) {}
  
  // arithmetic
  NumberFieldElement<R> operator-() const;
  NumberFieldElement<R> operator+(const NumberFieldElement<R> & ) const;
  NumberFieldElement<R> operator-(const NumberFieldElement<R> & ) const;
  NumberFieldElement<R> operator*(const NumberFieldElement<R> & ) const;
  NumberFieldElement<R> operator/(const NumberFieldElement<R> & ) const;

  NumberFieldElement<R> operator*(const R & ) const;
  NumberFieldElement<R> operator/(const R & ) const;

  NumberFieldElement<R> & operator+=(const NumberFieldElement<R> & );
  NumberFieldElement<R> & operator-=(const NumberFieldElement<R> & );
  NumberFieldElement<R> & operator*=(const NumberFieldElement<R> & );
  NumberFieldElement<R> & operator/=(const NumberFieldElement<R> & );

  NumberFieldElement<R>& operator*=(const R & );
  NumberFieldElement<R>& operator/=(const R & );


  static init_modulus(const Polynomial<R> & mod)
  {this->modulus = mod; }
  
protected:
  static UnivariatePoly<R> modulus;

  UnivariatePoly< Rational<R> > elt;
};

#include "NumberField.inl"

#endif // _NUMBER_FIELD_H
