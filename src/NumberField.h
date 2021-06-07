#ifndef __NUMBER_FIELD_H_
#define __NUMBER_FIELD_H_

#include "Polynomial.h"
#include "Rational.h"

template<typename R>
class NumberField :public std::enable_shared_from_this< const NumberField<R> > {
public:
  NumberField(const UnivariatePoly< Rational<R> > & mod) : f(mod) {}

  const UnivariatePoly< Rational<R> > & modulus() const
  {return this->f; }

  std::shared_ptr< const NumberField<R> > getptr() const
  {
    return std::enable_shared_from_this< const NumberField<R> >::shared_from_this();
  }
  
protected:
  UnivariatePoly< Rational<R> > f;
};

template<typename R>
class NumberFieldElement
{
public:

  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld) : K(fld) {}
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const UnivariatePoly< Rational<R> > & poly)
    : K(fld), elt(poly) {}
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const R & a)
    : K(fld), elt(a) {} 
  NumberFieldElement(std::shared_ptr<const NumberField<R> > fld,
		     const Rational<R> & a) : K(fld), elt(a) {}
  
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

  NumberFieldElement<R> inverse(void) const;
  
protected:
  std::shared_ptr<const NumberField<R> > K;
  
  UnivariatePoly< Rational<R> > elt;
};

#include "NumberField.inl"

#endif // _NUMBER_FIELD_H
