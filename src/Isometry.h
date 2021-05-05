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
  Isometry() : a(SquareMatrix<R, n>::identity()) {}

  Isometry(const SquareMatrix<R, n> & mat) : a(mat) {}

  // access - set/get
  void set_values(const SquareMatrix<R, n> & mat)
  { this->a = mat; }

  void set_identity(void)
  { this->a = SquareMatrix<R, n>::identity(); }

  const R & operator()(size_t i, size_t j) const
  { return this->a(i, j); }

  R & operator()(size_t i, size_t j)
  { return this->a(i, j); }

  // basic operations
  
  Isometry<R, n> inverse(void) const
  { return Isometry(this->a.inverse());}

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

  SquareMatrix<R, n> a;
};

#include "Isometry.inl"

#endif // __ISOMETRY_H_
