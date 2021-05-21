#ifndef __POLYNOMIAL_H_
#define __POLYNOMIAL_H_

#include "Fp.h"
#include "SquareMatrix.h"

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
  PolynomialFp(std::shared_ptr<const Fp<R,S>> GF, size_t i);

  // create a polynomial from a bilinear form
  template<size_t n>
  PolynomialFp(const SquareMatrixFp<R, S, n> & );
  
  // access

  // get methods
  const FpElement<R, S> & const_coefficient() const;
  // coefficient of x_i
  const FpElement<R, S> & coefficient(size_t i) const;
  // coefficient of x_i x_j
  const FpElement<R, S> & coefficient(size_t i, size_t j) const;

  const FpElement<R, S> &
  coefficient(const std::multiset<size_t> & mon) const;

  const std::map< std::multiset<size_t>, FpElement<R,S> > & monomials() const
  {return mons; }

  PolynomialFp<R,S> quadratic_part() const;
  std::vector< FpElement<R,S> > linear_part(size_t rank) const;

  // conversion, assignment operator
  PolynomialFp<R,S> & operator=(const PolynomialFp<R,S> & );
  PolynomialFp<R,S> & operator=(const FpElement<R,S> & );
  
  // arithmetic
  PolynomialFp<R,S> operator-() const;
  PolynomialFp<R,S> operator+(const PolynomialFp<R,S> & other) const;
  PolynomialFp<R,S> operator-(const PolynomialFp<R,S> & other) const;
  PolynomialFp<R,S> operator*(const PolynomialFp<R,S> & other) const;
  PolynomialFp<R,S> operator*(const FpElement<R,S> & a) const;

  FpElement<R,S> evaluate(const std::vector< FpElement<R,S> > & vec) const;
  PolynomialFp<R,S> evaluate(const std::vector<PolynomialFp<R,S> > & vec) const;

  // booleans
  bool operator==(const PolynomialFp<R, S> & other) const;
  bool operator!=(const PolynomialFp<R, S> & other) const;
  bool operator==(const FpElement<R, S> & other) const;
  bool operator!=(const FpElement<R, S> & other) const;

  
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
