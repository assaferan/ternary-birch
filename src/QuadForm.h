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
  QuadForm_Base() : is_reduced_(false) {}
  // from a vector of n(n+1)/2 elements
  QuadForm_Base(const SymVec& coeffs);
  QuadForm_Base(const SquareMatrix<R,n> & B) :
    B_(B), is_reduced_(false) {}
  
  // assignment
  QuadForm_Base<R,n>& operator=(const QuadForm_Base<R,n> &);

  // access
  R discriminant(void) const;

  bool operator==(const QuadForm_Base<R, n>& q) const
  { return (this->B_ == q.B_); }

  bool operator!=(const QuadForm_Base<R, n>& q) const
  {return !((*this)==q);}

  // !! TODO - check if there is a factor 2 here !!
  R evaluate(const Vector<R, n>& vec) const
  { return (vec, (this->B_) * vec) / 2; }

  const SquareMatrix<R, n> & bilinear_form() const
  { return this->B_; }

  bool is_reduced() const { return this->is_reduced_; }

  size_t num_automorphisms() const;
  
  Vector<R, n> orthogonalize_gram() const;

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
    SquareMatrix< FpElement<S, T>, n> q_mod;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	q_mod(i,j) = GF->mod(this->B_(i,j));
    QuadFormFp<S,T,n> q(q_mod);
    return q;
  }

  static QuadForm<R,n> reduce(const QuadForm<R,n> & q,
			      Isometry<R,n> & isom);

  std::set<Isometry<R, n>> proper_automorphisms() const;
  
protected:
  // a more general approach - the matrix representing the
  // bilinear form Q(x+y)-Q(x)-Q(y) (so Q(v) = 1/2 B(v,v))
  SquareMatrix<R, n> B_;
  bool is_reduced_;
  // we save it for quick access
  size_t num_aut_;

  // helper functions
  
  // reduce the form to a Minkowski reduced form
  // This is non-constant because we update the members
  // updates also the automorphism group of the lattice
 
  static size_t i_reduce(SquareMatrix<R, n> & qf,
			 Isometry<R,n> & isom,
			 std::set< Isometry<R, n> > & auts);

  static bool permutation_reduction(SquareMatrix<R, n> & qf,
				    Isometry<R,n> & isom,
				    std::set< Isometry<R, n> > & auts);
  
  static bool sign_normalization(SquareMatrix<R, n> & qf,
				 Isometry<R,n> & isom,
				 std::set< Isometry<R, n> > & auts);
  
  static bool norm_echelon(SquareMatrix<R, n> & qf, Isometry<R,n> & isom);
  
  static bool neighbor_reduction(SquareMatrix<R, n> & qf,
				 Isometry<R,n> & isom,
				 std::set< Isometry<R, n> > & auts);
  
  static size_t generate_auts(std::set< Isometry<R, n> > & auts);
  
  // static helper functions

  static std::vector< std::vector<size_t> > all_perms(size_t m);

  static void greedy(SquareMatrix<R,n>& q, Isometry<R,n>& s);
  
  static int hasse(const Vector<R, n>& , const R & );
  
  static Vector<R, n> voronoi_bounds(void);

  // update in-place q and iso according to the closest vector
  // to the space spanned by the n-1 first ones
  static void closest_lattice_vector(SquareMatrix<R,n> &q,
				     Isometry<R,n> & iso);
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
  QuadFormFp(const SquareMatrix< FpElement<R, S>, n> & mat) : 
    QuadForm< FpElement<R, S>, n>(mat) {this->GF = mat(0,0).field();}
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

  using QuadForm< FpElement<R, S> , n>::bilinear_form;
  using QuadForm< FpElement<R, S> , n>::discriminant;
  using QuadForm< FpElement<R, S> , n>::evaluate;

  R evaluate(const Vector<R, n>& vec) const {
    Vector< FpElement<R, S>, n> vec_mod = this->GF->mod(vec_mod);
    return (this->evaluate(vec_mod)).lift();
  }
  
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
