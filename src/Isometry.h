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

  Isometry(const SquareMatrix<R, n> & mat, const R & scale) :
    a(mat), scale(scale) {}

  // access - set/get
  const R & get_scale(void) const
  { return this->scale; }
  
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
  
  Isometry<R, n> inverse(void) const;

  Isometry<R, n> transpose(void) const
  { return Isometry(this->a.transpose(), this->scale); }

  // arithmetic
  
  Isometry<R, n> operator*(const Isometry<R, n>& s) const
  { return Isometry((this->a)*s.a, (this->scale)*s.scale); }

  Vector<R, n> operator*(const Vector<R, n>& vec) const
  { return (this->a)*vec; }

  // assignment
  Isometry<R, n> & operator=(const Isometry<R, n>& other)
  { if (this != &other) {
      this->a = other.a;
      this->scale = other.scale;
    }
    return *this;
  }
  
  SquareMatrix<R, n> transform(const SquareMatrix<R, n>& from) const;

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
