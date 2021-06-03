#ifndef __QUAD_FORM_H_
#define __QUAD_FORM_H_

#include <iostream>
#include <set>
#include <unordered_map>
#include <vector>

#include "birch.h"
#include "birch_util.h"
#include "Matrix.h"
#include "ParseNipp.h"
#include "Rational.h"

template<typename R, size_t n>
std::ostream& operator<<(std::ostream&, const QuadForm<R,n> &);

template<typename R, size_t n>
class QuadForm_Base
{
  public:
  typedef R SymVec[n*(n+1)/2];

  // c-tors
  QuadForm_Base() : is_reduced_(false), num_aut_init_(false) {}
  // from a vector of n(n+1)/2 elements
  QuadForm_Base(const SymVec& coeffs);
  QuadForm_Base(const SquareMatrix<R,n> & B) :
    B_(B), is_reduced_(false), num_aut_(0), num_aut_init_(false) {}
  
  // assignment
  QuadForm_Base<R,n>& operator=(const QuadForm_Base<R,n> &);

  // access
  R discriminant(void) const;

  bool operator==(const QuadForm_Base<R, n>& q) const
  { return (this->B_ == q.B_); }

  bool operator!=(const QuadForm_Base<R, n>& q) const
  {return !((*this)==q);}

  R evaluate(const Vector<R, n>& vec) const
  { return Vector<R, n>::inner_product(vec, (this->B_) * vec) / 2; }

  const SquareMatrix<R, n> & bilinear_form() const
  { return this->B_; }

  bool is_reduced() const { return this->is_reduced_; }

  size_t num_automorphisms() const;

  void set_num_aut(size_t num_aut)
  { this->num_aut_ = num_aut; this->num_aut_init_ = true; return;}

  void set_reduced()
  { this->is_reduced_ = true; return; }
  
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
  std::shared_ptr< QuadFormFp<S,T,n> > mod(std::shared_ptr< Fp<S,T> > GF) const
  {
    SquareMatrixFp<S, T, n> q_mod(GF);
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	q_mod(i,j) = GF->mod(this->B_(i,j));
    R p = GF->prime();
    if (p == 2) {
      for (size_t i = 0; i < n; i++) {
	R value = this->B_(i,i) / 2;
	q_mod(i,i) = GF->mod(value);
      }
    }
    std::shared_ptr< QuadFormFp<S,T,n> > q =
      std::make_shared< QuadFormFp<S, T, n> >(q_mod);
    return q;
  }

  static QuadForm<R,n> reduce(const QuadForm<R,n> & q,
			      Isometry<R,n> & isom,
			      bool calc_aut = false);

  std::set<Isometry<R, n>> proper_automorphisms() const;

  static std::vector< QuadForm<R, 5> > nipp_to_forms(NippEntry entry);
  
  // !! TODO - get_quinary_forms and nipp_to_forms should also work for
  // arbitrary R, no reason to restrict to Z, I think
  static std::vector<std::vector< QuadForm<R, 5> > >
  get_quinary_forms(const R & disc);

  static void greedy(SquareMatrix<R,n>& q, Isometry<R,n>& s, size_t dim = n);

  std::unordered_map< QuadForm<R, n>, Isometry<R,n> >
  generate_orbit() const;
  
protected:
  // a more general approach - the matrix representing the
  // bilinear form Q(x+y)-Q(x)-Q(y) (so Q(v) = 1/2 B(v,v))
  SquareMatrix<R, n> B_;
  bool is_reduced_;
  // we save it for quick access
  size_t num_aut_;
  bool num_aut_init_;

  // helper functions
  
  // reduce the form to a Minkowski reduced form
  // This is non-constant because we update the members
  // updates also the automorphism group of the lattice
 
  static size_t i_reduce(SquareMatrix<R, n> & qf,
			 Isometry<R,n> & isom,
			 std::set< Isometry<R, n> > & auts,
			 bool calc_aut = true);

  static bool permutation_reduction(SquareMatrix<R, n> & qf,
				    Isometry<R,n> & isom,
				    std::set< Isometry<R, n> > & auts,
				    bool calc_aut = true);
  
  static bool sign_normalization(SquareMatrix<R, n> & qf,
				 Isometry<R,n> & isom,
				 std::set< Isometry<R, n> > & auts,
				 bool calc_aut = true);
  
  static bool norm_echelon(SquareMatrix<R, n> & qf, Isometry<R,n> & isom);
  
  static bool neighbor_reduction(SquareMatrix<R, n> & qf,
				 Isometry<R,n> & isom,
				 std::set< Isometry<R, n> > & auts,
				 bool calc_aut = true);
  
  static size_t generate_auts(std::set< Isometry<R, n> > & auts);

  static Vector<R, n-1> voronoi_bounds(size_t dim = n);
  
  // static helper functions

  static std::vector< std::vector<size_t> > all_perms(size_t m);
  
  static int hasse(const Vector<R, n>& , const R & );

  // update in-place q and iso according to the closest vector
  // to the space spanned by the n-1 first ones
  static void closest_lattice_vector(SquareMatrix<R,n> &q,
				     Isometry<R,n> & iso,
				     size_t dim = n);

  static bool sign_normalization_slow(SquareMatrix<R, n> & qf,
				      Isometry<R,n> & isom,
				      std::set< Isometry<R, n> > & auts);

  static bool sign_normalization_fast(SquareMatrix<R, n> & qf,
				      Isometry<R,n> & isom);

  static std::vector<uint8_t> bit_transpose(const std::vector< uint8_t > & mat);
  static uint8_t bit_echelon_form(std::vector< uint8_t > & mat,
				  std::vector< uint8_t > & trans);
  
  static std::vector<uint8_t> kernel(const std::vector< uint8_t > & mat);

  std::unordered_map< QuadForm<R, n>, Isometry<R,n> >
  permutation_orbit() const;
  
  std::unordered_map<  QuadForm<R, n>, Isometry<R,n> >
  sign_orbit() const;
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

  //  friend std::ostream& operator<< <> (std::ostream&, const QuadForm&);

  using QuadForm_Base<R,n>::reduce;

  
protected:
 
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

  using QuadForm_Base<Z,n>::evaluate;
  using QuadForm_Base<Z,n>::reduce;
  
  static Z_QuadForm<3> get_quad_form(const std::vector<Z_PrimeSymbol>& input);

  using QuadForm_Base<Z,n>::generate_orbit;

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

  using QuadForm_Base<Z64,n>::evaluate;
  using QuadForm_Base<Z64,n>::reduce;
  using QuadForm_Base<Z64,n>::generate_orbit;
};

template<size_t n>
class QuadForm<Z128, n> : public QuadForm_Base<Z128,n>
{
public:
  QuadForm() : QuadForm_Base<Z128,n>() {}

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadForm(const typename QuadForm_Base<Z128,n>::SymVec& coeffs)
    : QuadForm_Base<Z128,n>(coeffs) {}

  QuadForm(const SquareMatrix<Z128, n> & B)
    : QuadForm_Base<Z128,n>(B) {}
 
  using QuadForm_Base<Z128,n>::operator==;
    
  using QuadForm_Base<Z128,n>::discriminant;
 
  W64 hash_value(void) const;

  using QuadForm_Base<Z128,n>::evaluate;
  using QuadForm_Base<Z128,n>::reduce;
  using QuadForm_Base<Z128,n>::generate_orbit;
};

template<typename R, typename S, size_t n>
class QuadFormFp : public QuadForm< FpElement<R, S> , n>
{
public:
  QuadFormFp(const SquareMatrix< FpElement<R, S>, n> & mat) : 
    QuadForm< FpElement<R, S>, n>(mat),
    GF(mat(0,0).field()),
    B_Fp(mat(0,0).field(), mat)
  {}

  QuadFormFp(const SquareMatrixFp<R, S, n> & mat) : 
    GF(mat.field()),
    B_Fp(mat.field(), mat)
  {}
  
  QuadFormFp(const typename QuadForm_Base<R,n>::SymVec& vec,
	     std::shared_ptr<const Fp<R,S>> GF) :
    QuadForm<FpElement<R, S> ,n>(GF->mod(vec)),
    GF(GF),
    B_Fp(GF, this->bilinear_form())
  {}

  const std::shared_ptr<const Fp<R,S>>& field(void) const
  {
    return this->GF;
  }

  const SquareMatrixFp<R, S, n> & bilinear_form() const
  { return B_Fp; }
  
  using QuadForm< FpElement<R, S> , n>::discriminant;
  
  FpElement<R, S> evaluate(const VectorFp<R, S, n>& v) const {
    R p = this->GF->prime();
    if (p == 2) return this->evaluate_p2(v);
    VectorFp<R, S, n> Bv = (this->bilinear_form()) * v;
    FpElement<R, S> two(GF, 2);
    return VectorFp<R, S, n>::inner_product(v, Bv)/two;
  }
  
  R evaluate(const Vector<R, n>& vec) const {
    VectorFp<R, S, n> v = this->GF->mod(vec);
    return (this->evaluate(v)).lift();
  }
  
  bool isotropic_vector(VectorFp<R, S, n> &,
			size_t start = 0,
			bool deterministic = false) const;

  void decompose(SquareMatrixFp<R, S, n> &,
		 SquareMatrixFp<R, S, n> &,
		 bool deterministic = false) const;

protected:
  std::shared_ptr<const Fp<R,S>> GF;
  SquareMatrixFp<R, S, n> B_Fp;

  // To avoid unnecessary computation, we encode each of the three 2-isotropic
  // vectors as a coordinate of the return vector. Special care must be taken
  // to obtain the actual isotropic vectors when needed.
 
  bool isotropic_vector_p2(VectorFp<R, S, n> &, size_t start = 0) const;

  FpElement<R, S> evaluate_p2(const VectorFp<R, S, n>& v) const;

  void hyperbolize_form(SquareMatrixFp<R, S, n> &,
			SquareMatrixFp<R, S, n> &,
			bool deterministic = false,
			size_t start = 0) const;
  
  void split_hyperbolic_plane(const VectorFp<R, S, n> &,
			      SquareMatrixFp<R, S, n> &,
			      SquareMatrixFp<R, S, n> &,
			      size_t start = 0) const;
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
