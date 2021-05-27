#ifndef __NEIGHBOR_MANAGER_H_
#define __NEIGHBOR_MANAGER_H_

#include "birch.h"
#include "Polynomial.h"
#include "QuadForm.h"
#include "Isometry.h"
#include "Fp.h"

template<typename R, typename S, typename T, size_t n>
class NeighborManager
{
public:
  NeighborManager(const QuadForm<T, n>& q, std::shared_ptr<Fp<R,S>> GF,
		  size_t k = 1);

  void next_isotropic_subspace(void);

  inline GenusRep<T, n> get_reduced_neighbor_rep() const;

  // to representative of the line
  Vector<R, n> transform_vector(const GenusRep<T, n>& dst, Vector<R, n> src);
  
  void get_next_neighbor(void);

  QuadForm<T, n> build_neighbor(Isometry<T, n>& ) const;

  const std::vector< Vector<T, n> > & get_isotropic_subspace() const
  {return this->X; }

protected:
  std::shared_ptr< Fp<R,S> > GF;
  QuadForm<T, n> q;
  T disc;
  SquareMatrix< FpElement<R, S> , n> b;
  SquareMatrix<T, n> quot_gram;
  std::shared_ptr< SquareMatrixFp<R, S, n> > p_std_gram;
  std::shared_ptr< SquareMatrixFp<R, S, n> > p_basis;
  std::shared_ptr< PolynomialFp<R, S> > p_q_std;
  // dimension of the radical
  size_t rad_dim;
  // dimension of the anisotropic subspace
  size_t aniso_dim;
  // the Witt index (number of hyperbolic planes)
  size_t witt_index;

  VectorFp< R, S, n> vec;
  std::vector< std::vector< size_t> > pivots;
  size_t pivot_ptr;
  size_t k; // dimension of the isotropic subspace
  size_t skew_dim;
  std::shared_ptr< MatrixFp<R,S> > p_skew;
  std::vector<size_t> free_vars;
  std::vector<FpElement<R,S> > params;
  std::shared_ptr<Matrix<PolynomialFp<R,S> > > p_isotropic_param;
  std::vector< VectorFp<R, S, n> > iso_subspace;
  std::vector< Vector<T, n> > X, Z, U;

  // The 2-isotropic vectors were stored in binary within each of the
  // coordinates of `vec` and so we use this function to unpack them into
  // actual 2-isotropic vectors.
  Vector<R, n> isotropic_vector_p2(R t) const;

  // get all possible pivots
  static std::vector< std::vector<size_t> >
  __pivots(size_t dim, size_t aniso, size_t k);

  void __initialize_pivots(void);

  SquareMatrix<T,n> __gram(const SquareMatrix<T,n> & B, bool quot = true) const;
  
  void lift_subspace();
  void update_skew_space();
  void update_skew_matrix(size_t &, size_t &);
};

#include "NeighborManager.inl"

#endif // __NEIGHBOR_MANAGER_H
