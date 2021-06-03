#ifndef __POLYNOMIAL_H_
#define __POLYNOMIAL_H_

#include "Fp.h"
#include "SquareMatrix.h"

template<typename R>
class UnivariatePoly
{
public:
  // create the zero polynomial
  UnivariatePoly() {}
  // create the constant polynomial
  UnivariatePoly(const R &);

  // create polynomial from coefficients
  UnivariatePoly(const std::vector<R> &);
  
  // create the polynomial x
  static UnivariatePoly<R> x();
  
  // access
  // get methods
  R const_coefficient() const {return this->coefficient(0); }
  
  // coefficient of x^i
  R coefficient(size_t i) const;

  const std::vector<R> & coefficients() const
  {return this->coeffs; }

  // if poly == 0, returns -1
  int degree() const {return this->coeffs.size()-1; }

  // conversion, assignment operator
  UnivariatePoly<R> & operator=(const UnivariatePoly<R> & );
  UnivariatePoly<R> & operator=(const R & );
  
  // arithmetic
  UnivariatePoly<R> operator-() const;
  UnivariatePoly<R> operator+(const UnivariatePoly<R> & ) const;
  UnivariatePoly<R> operator-(const UnivariatePoly<R> & ) const;
  UnivariatePoly<R> operator*(const UnivariatePoly<R> & ) const;
  UnivariatePoly<R> operator*(const R & ) const;

  UnivariatePoly<R> & operator+=(const UnivariatePoly<R> & );
  UnivariatePoly<R> & operator-=(const UnivariatePoly<R> & );
  UnivariatePoly<R> & operator*=(const UnivariatePoly<R> & );
  UnivariatePoly<R>& operator*=(const R & );

  UnivariatePoly<R> evaluate(const UnivariatePoly<R> &) const;
  R evaluate(const R &) const;

  // booleans
  bool is_zero() const {return this->coeffs.is_empty();}
  bool operator==(const UnivariatePoly<R> & ) const;
  bool operator!=(const UnivariatePoly<R> & ) const;
  bool operator==(const R & ) const;
  bool operator!=(const R & ) const;
  
  
protected:
  std::vector<R> coeffs;
};

template<typename R>
UnivariatePoly<R> operator*(const R & a,
			    const UnivariatePoly<R>  & poly)
{ return poly*a; }

template<typename R>
std::ostream& operator<<(std::ostream&, const UnivariatePoly<R> &);

template<typename R, typename S>
class PolynomialFp
{
public:

  // create the zero polynomial
  PolynomialFp(std::shared_ptr<const Fp<R,S>> GF);
  // create the constant polynomial
  PolynomialFp(const FpElement<R,S> & a);
  PolynomialFp(std::shared_ptr<const Fp<R,S>> GF, const R & a);
  // create the polynomial x_i
  // PolynomialFp(std::shared_ptr<const Fp<R,S>> GF, size_t i);
  static PolynomialFp<R,S> x(std::shared_ptr<const Fp<R,S>> GF, size_t i);

  // create a polynomial from a bilinear form
  template<size_t n>
  PolynomialFp(const SquareMatrixFp<R, S, n> & );
  
  // access

  // get methods
  FpElement<R, S> const_coefficient() const;
  // coefficient of x_i
  FpElement<R, S> coefficient(size_t i) const;
  // coefficient of x_i x_j
  FpElement<R, S> coefficient(size_t i, size_t j) const;

  FpElement<R, S>
  coefficient(const std::multiset<size_t> & mon) const;

  const std::map< std::multiset<size_t>, FpElement<R,S> > & monomials() const
  {return mons; }

  PolynomialFp<R,S> quadratic_part() const;
  std::vector< FpElement<R,S> > linear_part(size_t rank) const;

  int degree(size_t i) const;

  // conversion, assignment operator
  PolynomialFp<R,S> & operator=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator=(const FpElement<R,S> & );
  
  // arithmetic
  PolynomialFp<R,S> operator-() const;
  PolynomialFp<R,S> operator+(const PolynomialFp<R,S> & ) const;
  PolynomialFp<R,S> operator-(const PolynomialFp<R,S> & ) const;
  PolynomialFp<R,S> operator*(const PolynomialFp<R,S> & ) const;
  PolynomialFp<R,S> operator*(const FpElement<R,S> & ) const;

  PolynomialFp<R,S> & operator+=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator-=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator*=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator*=(const FpElement<R,S> & );

  FpElement<R,S> evaluate(const std::vector< FpElement<R,S> > & ) const;
  PolynomialFp<R,S> evaluate(const std::vector<PolynomialFp<R,S> > & vec) const;

  // booleans
  bool is_zero() const;
  bool operator==(const PolynomialFp<R, S> & ) const;
  bool operator!=(const PolynomialFp<R, S> & ) const;
  bool operator==(const FpElement<R, S> & ) const;
  bool operator!=(const FpElement<R, S> & ) const;

  
protected:
  std::shared_ptr<const Fp<R,S>> GF;
  
  std::map< std::multiset<size_t>, FpElement<R,S> > mons;
  
};

template<typename R, typename S>
PolynomialFp<R,S> operator*(const FpElement<R,S> & a,
			      const PolynomialFp<R,S>  & poly)
{ return poly*a; }

template<typename R, typename S>
std::ostream& operator<<(std::ostream&, const PolynomialFp<R,S>&);

#include "Polynomial.inl"

#endif // __POLYNOMIAL_H_
