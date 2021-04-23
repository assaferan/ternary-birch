#ifndef __QUAD_FORM_H_
#define __QUAD_FORM_H_

#include "birch.h"
#include "birch_util.h"
#include "ParseNipp.h"

template<typename R, size_t Rank>
class QuadForm
{
public:
  typedef R RMat[Rank][Rank];
  typedef R RVec[Rank*(Rank+1)/2];
  
  QuadForm() = default;

  // a more general constructor
  // We adhere to magma convention - giving the rows
  // up to the diagonal
  QuadForm(const RVec& coeffs);

  // These are only relevant for 3, do something about it later on
  QuadForm(const R& a, const R& b, const R& c,
	   const R& f, const R& g, const R& h)
  {
    this->a_ = a; this->b_ = b; this->c_ = c;
    this->f_ = f; this->g_ = g; this->h_ = h;
  }
  
  const R& a(void) const { return this->a_; }
  const R& b(void) const { return this->b_; }
  const R& c(void) const { return this->c_; }
  const R& f(void) const { return this->f_; }
  const R& g(void) const { return this->g_; }
  const R& h(void) const { return this->h_; }
  
  R discriminant(void) const;

  bool operator==(const QuadForm<R>& q) const
  {
    return this->a_ == q.a_ && this->b_ == q.b_ && this->c_ == q.c_ &&
      this->f_ == q.f_ && this->g_ == q.g_ && this->h_ == q.h_;
  }

  W64 hash_value(void) const;

  R evaluate(const R& x, const R& y, const R& z) const
  {
    return x * (this->a_ * x + this->g_ * z + this->h_ * y) +
      y * (this->b_ * y + this->f_ * z) + z * z * this->c_;
  }

  R evaluate(const Vector3<R>& vec) const
  {
    return this->evaluate(vec.x, vec.y, vec.z);
  }

  inline const RMat & getBilinearForm() const
  {
    return this->B_;
  }

  template<typename S, typename T>
  QuadFormFp<S,T> mod(std::shared_ptr<Fp<S,T>> GF) const
  {
    QuadFormFp<S,T> q(GF->mod(this->a_), GF->mod(this->b_), GF->mod(this->c_),
		      GF->mod(this->f_), GF->mod(this->g_), GF->mod(this->h_), GF);
    return q;
  }

  static Z_QuadForm get_quad_form(const std::vector<PrimeSymbol<R>>& primes)
  {
    static_assert( std::is_same<R,Z>::value, "Implemented only for arbitrary precision types." );
    return Z_QuadForm(); // Make the compiler happy.
  }

  static std::vector< QuadForm<Z,5> >
  nippToForms(NippEntry);
  
  static std::vector< std::vector<QuadForm<Z,5> > >
  get_quinary_forms(const Z &);

  static int border(const QuadForm<R>&, int);

  static int num_automorphisms(const QuadForm<R>&);

  static const std::vector<Isometry<R>>&
  proper_automorphisms(const QuadForm<R>&);

  static QuadForm<R> reduce(const QuadForm<R>&, Isometry<R>&);

  friend std::ostream& operator<<(std::ostream&, const QuadForm&);

protected:
  // a more general approach - the matrix representing the
  // bilinear form Q(x+y)-Q(x)-Q(y) (so Q(v) = 1/2 B(v,v))
  RMat B_;
    
  R a_, b_, c_, f_, g_, h_;
};

template<typename R, typename S>
class QuadFormFp : public QuadForm<R>
{
public:
  QuadFormFp(const R& a, const R& b, const R& c,
	     const R& f, const R& g, const R& h,
	     std::shared_ptr<Fp<R,S>> GF) :
    QuadForm<R>(GF->mod(a), GF->mod(b), GF->mod(c),
		GF->mod(f), GF->mod(g), GF->mod(h))
  {
    this->GF = GF;
  }

  const std::shared_ptr<Fp<R,S>>& field(void) const
  {
    return this->GF;
  }

  R discriminant(void) const;
  R evaluate(const R& x, const R& y, const R& z) const;

  R evaluate(const Vector3<R>& vec) const
  {
    return this->evaluate(vec.x, vec.y, vec.z);
  }

  Vector3<R> isotropic_vector(void) const;

private:
  std::shared_ptr<Fp<R,S>> GF;

  // To avoid unnecessary computation, we encode each of the three 2-isotropic
  // vectors as a coordinate of the return vector. Special care must be taken
  // to obtain the actual isotropic vectors when needed.
  Vector3<R> isotropic_vector_p2(void) const;
};

namespace std
{
  template<typename R>
  struct hash<QuadForm<R>>
  {
    Z64 operator()(const QuadForm<R>& q) const
    {
      return q.hash_value();
    }
  };

  template<typename R>
  struct hash<GenusRep<R>>
  {
    Z64 operator()(const GenusRep<R>& rep) const
    {
      return rep.q.hash_value();
    }
  };
}

template<typename R>
std::ostream& operator<<(std::ostream& os, const Vector3<R>& vec)
{
  os << "Vector(" << vec.x << "," << vec.y << "," << vec.z << ")";
  return os;
}

template<typename R>
bool operator==(const Vector3<R>& vec1, const Vector3<R>& vec2)
{
  return vec1.x == vec2.x && vec1.y == vec2.y && vec1.z == vec2.z;
}

template<typename R>
Vector3<R> operator+(const Vector3<R>& a, const Vector3<R>& b)
{
  Vector3<R> res;
  res.x = a.x + b.x;
  res.y = a.y + b.y;
  res.z = a.z + b.z;
  return res;
}

template<>
Z_QuadForm Z_QuadForm::get_quad_form(const std::vector<Z_PrimeSymbol>& primes);

#include "QuadForm.inl"

#endif // __QUAD_FORM_H_
