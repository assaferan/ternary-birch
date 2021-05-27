#ifndef __ISOMETRY_H_
#define __ISOMETRY_H_

#include "birch.h"
#include "SquareMatrix.h"

// !! TODO - think if we really need a different class for Isometry

template<typename R, size_t n>
class Isometry
{
public:
  // c-tors
  Isometry() : a(SquareMatrix<R, n>::identity()), scale(Math<R>::one()) {}

  Isometry(const SquareMatrix<R, n> & mat) : a(mat), scale(Math<R>::one()) {}

  // access - set/get
  void set_values(const SquareMatrix<R, n> & mat)
  { this->a = mat; }

  void set_identity(void)
  { this->a.set_identity(); this->scale = Math<R>::one(); }

  void set_scale(const R & scale)
  { this->scale = scale; }

  const R & operator()(size_t i, size_t j) const
  { return this->a(i, j); }

  R & operator()(size_t i, size_t j)
  { return this->a(i, j); }

  // basic operations
  
  Isometry<R, n> inverse(void) const
  {
    // !! TODO - should be able to invert without using rationals
    // for example, can always track back (save the inverse for the ride)
    SquareMatrix< Rational<R>, n> a_rat;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++) {
	a_rat(i,j) = this->a(i,j);
	a_rat(i,j) /= this->scale;
      }
    a_rat = a_rat.inverse();
    SquareMatrix<R, n> a_inv;
    // Since this is an isometry, the inverse should be integral
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	a_inv(i,j) = (this->scale * a_rat(i,j)).floor();
#ifdef DEBUG
    R scale2 = (this->scale)*(this->scale);
    assert((a_inv * (this->a) == scale2*SquareMatrix<R, n>::identity()));
#endif
    return Isometry(a_inv);
  }

  Isometry<R, n> transpose(void) const
  { return Isometry(this->a.transpose()); }

  // arithmetic
  
  Isometry<R, n> operator*(const Isometry<R, n>& s) const
  { return Isometry((this->a)*s.a); }

  Vector<R, n> operator*(const Vector<R, n>& vec) const
  { return (this->a)*vec; }

  // assignment
  Isometry<R, n> & operator=(const Isometry<R, n>& other)
  { if (this != &other) {
      this->a = other.a;
    }
    return *this;
  }
  
  SquareMatrix<R, n> transform(const SquareMatrix<R, n>& from, R scalar) const;

  // we save some clocks by returning once a single coordinate is mismatched.
  bool is_isometry(const QuadForm<R, n>& from, const QuadForm<R, n>& to,
		   R scalar) const;
  
  void update_perm(const Vector<size_t, n> & perm);

  friend std::ostream& operator<<(std::ostream& os, const Isometry<R, n>& s)
  { os << s.a; return os; }

  bool operator==(const Isometry<R, n> & other) const
  {return (this->a == other.a);}

  bool operator<(const Isometry<R, n> & other) const
  {return (this->a < other.a);}

  SquareMatrix<R, n> a;
  
protected:
 
  R scale;
};

#include "Isometry.inl"

#endif // __ISOMETRY_H_
