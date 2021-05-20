#ifndef __NEIGHBOR_MANAGER_H_
#define __NEIGHBOR_MANAGER_H_

#include "birch.h"
#include "QuadForm.h"
#include "Isometry.h"
#include "Fp.h"

template<typename R, typename S, typename T, size_t n>
class NeighborManager
{
public:
  NeighborManager(const QuadForm<T, n>& q, std::shared_ptr<Fp<R,S>> GF, size_t k = 1);

  Vector<R, n> isotropic_vector(R t) const;

  inline GenusRep<T, n> get_reduced_neighbor_rep(R t) const;

  // to representative of the line
  Vector<R, n> transform_vector(const GenusRep<T, n>& dst, Vector<R, n> src);
  
  QuadForm<T, n> get_neighbor(R t, Isometry<T, n>& s) const;

  QuadForm<T, n> build_neighbor(Vector<R, n>& vec2, Isometry<T, n>& s) const;

protected:
  std::shared_ptr<Fp<R,S>> GF;
  QuadForm<T, n> q;
  T disc;
  SquareMatrix< FpElement<R, S> , n> b;
  std::shared_ptr< SquareMatrixFp<R, S, n> > p_std_gram;
  std::shared_ptr< SquareMatrixFp<R, S, n> > p_basis;
  // dimension of the radical
  size_t rad_dim;
  // dimension of the anisotropic subspace
  size_t aniso_dim;
  // the Witt index (number of hyperbolic planes)
  size_t witt_index;

  VectorFp< R, S, n> vec;
  std::vector< std::vector< size_t> > pivots;

  // The 2-isotropic vectors were stored in binary within each of the
  // coordinates of `vec` and so we use this function to unpack them into
  // actual 2-isotropic vectors.
  Vector<R, n> isotropic_vector_p2(R t) const;

  // get all possible pivots
  static std::vector< std::vector<size_t> >
  __pivots(size_t dim, size_t aniso, size_t k);
};

#include "NeighborManager.inl"

#endif // __NEIGHBOR_MANAGER_H
