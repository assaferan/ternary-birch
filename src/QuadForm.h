#ifndef __QUAD_FORM_H_
#define __QUAD_FORM_H_

#include <iostream>
#include <set>
#include <vector>

#include "birch.h"
#include "birch_util.h"
#include "Matrix.h"
#include "ParseNipp.h"
#include "Rational.h"

template<typename R, size_t n>
std::ostream& operator<<(std::ostream&, const QuadForm<R,n>&);

template<typename R, size_t n>
class QuadForm_Base
{
  public:
  typedef R SymVec[n*(n+1)/2];

  // c-tors
  QuadForm_Base() : this->is_aut_init_(false) {}
  // from a vector of n(n+1)/2 elements
  QuadForm_Base(const SymVec& coeffs);
  QuadForm_Base(const SquareMatrix<R,n> & B) :
    this->B_(B), this->is_aut_init_(false) {}
  
  // assignment
  QuadForm_Base<R,n>& operator=(const QuadForm_Base<R,n> & other)
  {
    if ((*this) != other) {
      this->B_ = other.B_;
      this->is_aut_init_ = other.is_aut_init_;
      if (is_aut_init_) {
	this->aut_ = other.aut_;
	this->B_red_ = other.B_red_;
	this->isom_ = other.isom_;
      }
    }
    return *this;
  }

  // access
  R discriminant(void) const;

  bool operator==(const QuadForm_Base<R, n>& q) const
  { return (this->B_ == q.B_); }

  bool operator!=(const QuadForm_Base<R, n>& q) const
  {return !((*this)==q);}

  // !! TODO - check if there is a factor 2 here !!
  R evaluate(const Vector<R, n>& vec) const
  { return (vec, (this->B_) * vec) / 2; }

  inline const RMat & bilinear_form() const
  { return this->B_; }

  std::vector<R> orthogonalize_gram() const;

  R invariants(std::set<R> & , size_t& ) const;
  
  R invariants(std::set<std::pair<R, int> > &, size_t& ) const;

  struct jordan_data {
    std::vector< Matrix< Rational<R> > > matrices;
    std::vector< Matrix< Rational<R> > > grams;
    std::vector<size_t> exponents;
  };
  
  jordan_data jordan_decomposition(const R & p) const;

  template<typename S, typename T>
  QuadFormFp<S,T,n> mod(std::shared_ptr<Fp<S,T>> GF) const
  {
    QuadFormFp<S,T,n> q(GF->mod(this->a_), GF->mod(this->b_), GF->mod(this->c_),
		      GF->mod(this->f_), GF->mod(this->g_), GF->mod(this->h_), GF);
    return q;
  }

  size_t num_automorphisms()
  {if (!is_aut_init_) reduce(); return aut_.size(); }

  // reduce the form to a Minkowski reduced form
  // This is non-constant because we update the members 
  void reduce(void);
  
protected:
  // a more general approach - the matrix representing the
  // bilinear form Q(x+y)-Q(x)-Q(y) (so Q(v) = 1/2 B(v,v))
  SquareMatrix<R, n> B_;
  SquareMatrix<R, n> B_red_;
  Isometry<R, n> isom_;
  std::set< Isometry<R,n> > aut_;
  bool is_aut_init_;

  // helper functions
  static int Hasse(const std::vector<R>& , const R & );
  
};

template<typename R, size_t n>
class QuadForm : public QuadForm_Base<R, n>
{
public:
  QuadForm() : QuadForm_Base<R,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadForm(const typename QuadForm_Base<R,n>::SymVec& coeffs)
    : QuadForm_Base<R,n>(coeffs) {}

  QuadForm(const SquareMatrix<R, n> & B)
    : QuadForm_Base<R,n>(B) {}

  friend std::ostream& operator<< <> (std::ostream&, const QuadForm&);

  using QuadForm_Base<R,n>::reduce;
  
};

template<size_t n>
class QuadForm<Z, n> : public QuadForm_Base<Z,n>
{
public:
  QuadForm() : QuadForm_Base<Z,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadForm(const typename QuadForm_Base<Z,n>::SymVec& coeffs)
    : QuadForm_Base<Z,n>(coeffs) {}

  QuadForm(const SquareMatrix<Z, n> & B)
    : QuadForm_Base<Z,n>(B) {}

  using QuadForm_Base<Z,n>::operator==;
    
  using QuadForm_Base<Z,n>::discriminant;
  
  W64 hash_value(void) const;

  using QuadForm_Base<Z, n>::evaluate;
  using QuadForm_Base<Z,n>::reduce;

  // !! TODO - get_quinary_forms and nipp_to_forms should also work for
  // arbitrary R, no reason to restrict to Z, I think
  static std::vector<std::vector< Z_QuadForm<5> > >
  get_quinary_forms(const Z & disc);

  static Z_QuadForm<3> get_quad_form(const std::vector<Z_PrimeSymbol>& input);

  static std::vector< Z_QuadForm<5> > nipp_to_forms(NippEntry entry);
};

// Check which ones I really need
template<size_t n>
class QuadForm<Z64, n> : public QuadForm_Base<Z64,n>
{
public:
  QuadForm() : QuadForm_Base<Z64,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadForm(const typename QuadForm_Base<Z64,n>::SymVec& coeffs)
    : QuadForm_Base<Z64,n>(coeffs) {}

  QuadForm(const SquareMatrix<Z64, n> & B)
    : QuadForm_Base<Z64,n>(B) {}
 
  using QuadForm_Base<Z64,n>::operator==;
    
  using QuadForm_Base<Z64,n>::discriminant;
 
  W64 hash_value(void) const;

  using QuadForm_Base<Z64, n>::evaluate;
  using QuadForm_Base<Z64, n>::reduce;
};

template<typename R, typename S, size_t n>
class QuadFormFp : public QuadForm< FpElement<R, S> , n>
{
public:
  QuadFormFp(const typename QuadForm_Base<R,n>::SymVec& vec,
	     std::shared_ptr<Fp<R,S>> GF) :
    QuadForm<FpElement<R, S> ,n>(GF->mod(vec))
  {
    this->GF = GF;
  }

  const std::shared_ptr<Fp<R,S>>& field(void) const
  {
    return this->GF;
  }

  using QuadForm< FpElement<R, S> , n>::discriminant;
  using QuadForm< FpElement<R, S> , n>::evaluate;

  Vector<FpElement<R,S>,n> isotropic_vector(void) const;

protected:
  std::shared_ptr<Fp<R,S>> GF;

  // To avoid unnecessary computation, we encode each of the three 2-isotropic
  // vectors as a coordinate of the return vector. Special care must be taken
  // to obtain the actual isotropic vectors when needed.
 
  Vector<FpElement<R,S>, n> isotropic_vector_p2(void) const;
};

namespace std
{
  template<typename R, size_t n>
  struct hash<QuadForm<R,n>>
  {
    Z64 operator()(const QuadForm<R,n>& q) const
    {
      return q.hash_value();
    }
  };

  template<typename R, size_t n>
  struct hash<GenusRep<R, n>>
  {
    Z64 operator()(const GenusRep<R, n>& rep) const
    {
      return rep.q.hash_value();
    }
  };
}

// for some reason can't override operator<< here
template <typename R>
void pretty_print(std::ostream & os,std::vector<R> vec)
{
  for (size_t i = 0; i < vec.size(); i++)
    os << vec[i] << " ";
  os << std::endl;
  return;
}

#include "QuadForm.inl"

#endif // __QUAD_FORM_H_
