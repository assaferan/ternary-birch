#ifndef __POLYNOMIAL_H_
#define __POLYNOMIAL_H_

#include "Fp.h"
#include "SquareMatrix.h"

#include <set>
#include <unordered_map>

template<typename R, typename S>
class UnivariatePolyFp;

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

  // create the polynomial x^i
  static UnivariatePoly<R> x(size_t i = 1);
  
  // access
  // get methods
  R const_coefficient() const {return this->coefficient(0); }
  
  // coefficient of x^i
  R coefficient(size_t i) const;

  // leading coefficient
  R lead() const { return this->coefficient(this->degree()); }

  const std::vector<R> & coefficients() const
  {return this->coeffs; }

  // if poly == 0, returns -1
  int degree() const {return this->coeffs.size()-1; }

  R content() const;

  // conversion, assignment operator
  template<typename T>
  UnivariatePoly<R> & operator=(const UnivariatePoly<T> & );
  UnivariatePoly<R> & operator=(const R & );
  
  // arithmetic
  UnivariatePoly<R> operator-() const;
  UnivariatePoly<R> operator+(const UnivariatePoly<R> & ) const;
  UnivariatePoly<R> operator-(const UnivariatePoly<R> & ) const;
  UnivariatePoly<R> operator*(const UnivariatePoly<R> & ) const;
  UnivariatePoly<R> operator/(const UnivariatePoly<R> & ) const;
  UnivariatePoly<R> operator%(const UnivariatePoly<R> & ) const;
  UnivariatePoly<R> operator*(const R & ) const;
  UnivariatePoly<R> operator/(const R & ) const;
  UnivariatePoly<R> operator%(const R & ) const;

  UnivariatePoly<R> & operator+=(const UnivariatePoly<R> & );
  UnivariatePoly<R> & operator-=(const UnivariatePoly<R> & );
  UnivariatePoly<R> & operator*=(const UnivariatePoly<R> & );
  UnivariatePoly<R> & operator/=(const UnivariatePoly<R> & );
  UnivariatePoly<R> & operator%=(const UnivariatePoly<R> & );
  UnivariatePoly<R>& operator*=(const R & );
  UnivariatePoly<R>& operator/=(const R & );
  UnivariatePoly<R>& operator%=(const R & );

  template<typename S, typename T>
  UnivariatePolyFp<S, T> mod(std::shared_ptr< const Fp<S, T> >) const;
  
  UnivariatePoly<R> evaluate(const UnivariatePoly<R> &) const;
  R evaluate(const R &) const;
  template<typename S>
  Matrix<S> evaluate(const Matrix<S> &) const;

  // booleans
  bool is_zero() const {return this->coeffs.empty();}
  bool operator==(const UnivariatePoly<R> & ) const;
  bool operator!=(const UnivariatePoly<R> & ) const;
  bool operator==(const R & ) const;
  bool operator!=(const R & ) const;

  // hash
  W64 hash_value(void) const;
  
  // algorithms
  UnivariatePoly<R> derivative() const;

  static void div_rem(const UnivariatePoly<R> & f,
		      const UnivariatePoly<R> & g,
		      UnivariatePoly<R> & q,
		      UnivariatePoly<R> & r);
  
  static UnivariatePoly<R> gcd(const UnivariatePoly<R> &,
			       const UnivariatePoly<R> &);

  static UnivariatePoly<R> xgcd(const UnivariatePoly<R> & f,
				const UnivariatePoly<R> & g,
				UnivariatePoly<R> & s,
				UnivariatePoly<R> & t);
  
  std::unordered_map< UnivariatePoly<R>, size_t > factor() const;
  
protected:
  std::vector<R> coeffs;

  void eliminate_deg();
  
  // these helper methods are needed for factorization
  
  template<typename S, typename T>
  void hensel_step(std::vector<UnivariatePoly<R> > & u,
		   std::vector<UnivariatePoly<R> > & v,
		   std::shared_ptr< const Fp<S,T> > GF,
		   size_t i) const;

  template<typename S, typename T>
  std::vector< UnivariatePoly<R> >
  hensel_lift(const std::vector<UnivariatePolyFp<S, T> > & g,
	      size_t a) const;

  std::vector< UnivariatePoly<R> >
  trial_factor(const std::vector<UnivariatePoly< R > > & u,
	       const R & N) const;

  void
  find_trial_factor(const std::vector< UnivariatePoly<R> > & u,
		    const R & N,
		    size_t & j,
		    std::set<size_t> & C,
		    size_t & s,
		    std::vector< UnivariatePoly<R> > & gs);

  R landau_mignotte() const;

  std::vector< UnivariatePoly<R> > squarefree_factor() const;

  static std::set< std::set<size_t> >
  subsets(const std::set<size_t> & S, size_t k);
};

template<typename R>
UnivariatePoly<R> operator*(const R & a,
			    const UnivariatePoly<R>  & poly)
{ return poly*a; }

template<typename R>
std::ostream& operator<<(std::ostream&, const UnivariatePoly<R> &);

namespace std
{
  template<typename R>
  struct hash< UnivariatePoly<R> >
  {
    Z64 operator()(const UnivariatePoly<R>& p) const
    {
      return p.hash_value();
    }
  };
}

template<typename R, typename S>
class UnivariatePolyFp : public UnivariatePoly< FpElement<R,S> >
{
public:
  UnivariatePolyFp(std::shared_ptr< const Fp<R,S>> GF)
  { this->GF_ = GF; }

  // create the constant polynomial
  UnivariatePolyFp(const FpElement<R, S> & a);
  
  // create polynomial from coefficients
  UnivariatePolyFp(const std::vector< FpElement<R,S> > & v)
    : UnivariatePoly< FpElement<R,S> >(v)
  {this->GF_ = v[0].field(); }

  UnivariatePolyFp(const UnivariatePoly< FpElement<R, S> > & other);
  
  // create the polynomial x^i
  static UnivariatePolyFp<R,S> x(std::shared_ptr< const Fp<R,S>> GF,
				 size_t i = 1);
  
  // access
  const std::shared_ptr< const Fp<R,S> > & field() const
  {return this->GF_;}

  FpElement<R,S> coefficient(size_t i) const;
  FpElement<R,S> content() const;

  // arithmetic
  UnivariatePolyFp<R, S> operator*(const UnivariatePolyFp<R,S> & ) const;
  UnivariatePolyFp<R, S> operator/(const UnivariatePolyFp<R,S> & ) const;
  UnivariatePolyFp<R, S> operator%(const UnivariatePolyFp<R,S> & ) const;

  UnivariatePolyFp<R,S> & operator*=(const UnivariatePolyFp<R,S> & );
  UnivariatePolyFp<R,S> & operator/=(const UnivariatePolyFp<R,S> & );
  UnivariatePolyFp<R,S> & operator%=(const UnivariatePolyFp<R,S> & );

  UnivariatePolyFp<R,S> derivative() const;
  
  std::vector< UnivariatePolyFp<R,S> > sqf_factor() const;

  UnivariatePoly<R> lift() const;
  
  UnivariatePolyFp<R,S>
  pow_mod(size_t, const UnivariatePolyFp<R,S> & ) const;

  // assignment and conversion
  using UnivariatePoly< FpElement<R,S> >::operator=;
  UnivariatePolyFp<R,S> & operator=(const UnivariatePoly< FpElement<R,S> > &);
  
  
  // couldn't make it work with inheritance
  static UnivariatePolyFp<R,S> gcd(const UnivariatePolyFp<R,S> & f,
				   const UnivariatePolyFp<R,S> & g);
				    
  static UnivariatePolyFp<R,S> xgcd(const UnivariatePolyFp<R,S> & f,
				    const UnivariatePolyFp<R,S> & g,
				    UnivariatePolyFp<R,S> & s,
				    UnivariatePolyFp<R,S> & t);

  static void div_rem(const UnivariatePolyFp<R,S> & f,
		      const UnivariatePolyFp<R,S> & g,
		      UnivariatePolyFp<R,S> & q,
		      UnivariatePolyFp<R,S> & r);
  
  /*
  using UnivariatePoly< FpElement<R,S> >::div_rem;
  using UnivariatePoly< FpElement<R,S> >::gcd;
  using UnivariatePoly< FpElement<R,S> >::xgcd;
  */
  
protected:
  std::shared_ptr< const Fp<R,S>> GF_;
  
  std::vector< UnivariatePolyFp<R,S> >
  cz_eq_deg_partial_factor(size_t r) const;

  std::vector< UnivariatePolyFp<R,S> > cz_eq_deg_factor(size_t r) const;

  std::vector< UnivariatePolyFp<R,S> > cz_distinct_deg_factor() const;
  
};

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
