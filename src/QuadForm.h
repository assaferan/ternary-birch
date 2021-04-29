#ifndef __QUAD_FORM_H_
#define __QUAD_FORM_H_

#include <iostream>
#include <set>
#include <vector>

#include "birch.h"
#include "birch_util.h"
#include "Matrix.h"
#include "ParseNipp.h"

template<typename R, size_t n>
std::ostream& operator<<(std::ostream&, const QuadForm<R,n>&);

template<typename R, size_t n>
class QuadForm
{
public:
  typedef R RMat[n][n];
  typedef R RVec[n*(n+1)/2];
  typedef R RDiag[n];
  
  QuadForm() = default;

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadForm(const RVec& coeffs);

  QuadForm(const RMat& B)
  {
    for (size_t row = 0; row < n; row++)
      for (size_t col = 0; col < n; col++)
	this->B_[row][col] = B[row][col];
  }

  // These are only relevant for 3, do something about it later on
  QuadForm(const R& a, const R& b, const R& c,
	   const R& f, const R& g, const R& h)
  {
    this->a_ = a; this->b_ = b; this->c_ = c;
    this->f_ = f; this->g_ = g; this->h_ = h;
  }

  // assignment
  QuadForm<R,n>& operator=(const QuadForm<R,n> & other)
  {
    if ((*this) != other) {
      for (size_t i = 0; i < n; i++)
	for (size_t j = 0; j < n; j++)
	  this->B_[i][j] = other.B_[i][j];
    }
    return *this;
  }
  
  const R& a(void) const { return this->a_; }
  const R& b(void) const { return this->b_; }
  const R& c(void) const { return this->c_; }
  const R& f(void) const { return this->f_; }
  const R& g(void) const { return this->g_; }
  const R& h(void) const { return this->h_; }

  // access
  R discriminant(void) const;

  bool operator==(const QuadForm<R, n>& q) const
  {
    return this->a_ == q.a_ && this->b_ == q.b_ && this->c_ == q.c_ &&
      this->f_ == q.f_ && this->g_ == q.g_ && this->h_ == q.h_;
  }

  bool operator!=(const QuadForm<R, n>& q) const
  {return !((*this)==q);}

  W64 hash_value(void) const;

  R evaluate(const R& x, const R& y, const R& z) const
  {
    return x * (this->a_ * x + this->g_ * z + this->h_ * y) +
      y * (this->b_ * y + this->f_ * z) + z * z * this->c_;
  }

  R evaluate(const Vector<R, n>& vec) const
  {
    return this->evaluate(vec.x, vec.y, vec.z);
  }

  inline const RMat & bilinear_form() const
  {
    return this->B_;
  }

  std::vector<R> orthogonalize_gram() const;

  R invariants(std::set<R> & , size_t& ) const;
  
  R invariants(std::set<std::pair<R, int> > &, size_t& ) const;

  struct jordan_data {
    std::vector< Matrix<R> > matrices;
    std::vector< Matrix<R> > grams;
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

  static int border(const QuadForm<R, n>&, int);

  static int num_automorphisms(const QuadForm<R, n>&);

  static const std::vector<Isometry<R,n>>&
  proper_automorphisms(const QuadForm<R, n>&);

  static QuadForm<R, n> reduce(const QuadForm<R, n>&, Isometry<R,n>&);

  friend std::ostream& operator<< <> (std::ostream&, const QuadForm&);

protected:
  // a more general approach - the matrix representing the
  // bilinear form Q(x+y)-Q(x)-Q(y) (so Q(v) = 1/2 B(v,v))
  RMat B_;
    
  R a_, b_, c_, f_, g_, h_;

private:
  
  static int Hasse(const std::vector<R>& , const R & );
  static R inner_product(const RMat & F, const RMat & S,
		  size_t idx1, size_t idx2);
};

template<typename R, typename S, size_t n>
class QuadFormFp : public QuadForm<R, n>
{
public:
  QuadFormFp(const R& a, const R& b, const R& c,
	     const R& f, const R& g, const R& h,
	     std::shared_ptr<Fp<R,S>> GF) :
    QuadForm<R,n>(GF->mod(a), GF->mod(b), GF->mod(c),
		GF->mod(f), GF->mod(g), GF->mod(h))
  {
    this->GF = GF;
  }

  const std::shared_ptr<Fp<R,S>>& field(void) const
  {
    return this->GF;
  }

  R discriminant(void) const;

  R evaluate(const Vector<R, n>& vec) const
  {
    // stub - !!! TODO !! complete
    return 0;
  }

  Vector<R,n> isotropic_vector(void) const;

protected:
  std::shared_ptr<Fp<R,S>> GF;

  // To avoid unnecessary computation, we encode each of the three 2-isotropic
  // vectors as a coordinate of the return vector. Special care must be taken
  // to obtain the actual isotropic vectors when needed.
  Vector3<R> isotropic_vector_p2(void) const;
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

template<typename R, size_t n>
std::ostream& operator<<(std::ostream& os, const Vector<R,n>& vec)
{
  os << "Vector(" << vec.x << "," << vec.y << "," << vec.z << ")";
  return os;
}

template<typename R, size_t n>
bool operator==(const Vector<R,n>& vec1, const Vector<R,n>& vec2)
{
  return vec1.x == vec2.x && vec1.y == vec2.y && vec1.z == vec2.z;
}

template<typename R, size_t n>
Vector<R, n> operator+(const Vector<R, n>& a, const Vector<R,n>& b)
{
  Vector<R,n> res;
  res.x = a.x + b.x;
  res.y = a.y + b.y;
  res.z = a.z + b.z;
  return res;
}

#include "QuadForm.inl"

#endif // __QUAD_FORM_H_
