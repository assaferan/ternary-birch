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
  NeighborManager(const QuadForm<T, n>& q, std::shared_ptr<Fp<R,S>> GF)
    : vec(GF)
  {
    this->q = q;
    this->disc = q.discriminant();

    std::shared_ptr<QuadFormFp<R,S,n> > qp = q.mod(GF);

    this->b = qp->bilinear_form();
	
    this->GF = GF;
    assert(qp->isotropic_vector(this->vec));

#ifdef DEBUG
    R prime = GF->prime();
    if (prime != 2) assert( qp.evaluate(vec) == 0 );
#endif

    SquareMatrixFp<R, S, n> gram(GF);
    SquareMatrixFp<R, S, n> basis(GF);
    qp->decompose(gram, basis);

#ifdef DEBUG
    std::cerr << "Performed Witt Decomposition on " << qp << ". ";
    std::cerr << "Resulting gram matrix is " << gram << ", ";
    std::cerr << "Resulting basis is " << basis << std::endl;
#endif
    
  }

  Vector<R, n> isotropic_vector(R t) const
  {
    Vector<R, n> res;

    R p = GF->prime();

    if (p == 2) return this->isotropic_vector_p2(t);

    // Stub
    // !! TODO - add code that generates the isotropic vector
    // corresponding to the parameter t

#ifdef DEBUG
    QuadFormFp<R,S,n> qp = this->q.mod(GF);
    assert( qp.evaluate(res) % this->GF->prime() == 0 );
#endif

    return res;
  }

  inline GenusRep<T, n> get_reduced_neighbor_rep(R t) const
  {
    GenusRep<T, n> rep;
    rep.q = this->get_neighbor(t, rep.s);
    rep.q = QuadForm<T, n>::reduce(rep.q, rep.s);
    return rep;
  }

  // to representative of the line
  Vector<R, n> transform_vector(const GenusRep<T, n>& dst, Vector<R, n> src)
  {
    Vector<T, n> temp;
    VectorFp<R, S, n> temp_mod = GF->mod(src);
    for (size_t i = 0; i < n; i++)
      temp[i] = temp_mod[i].lift();

    R p = GF->prime();

    // should that be inverse mod p?
    Isometry<T, n> sinv = dst.s.inverse();
    temp = sinv * temp;

#ifdef DEBUG
    for (size_t i = 0; i < n; i++)
      assert( temp[i] % p == 0 );
#endif

    for (size_t i = 0; i < n; i++)
      temp[i] /= p;

    Vector<R, n> vec;
    for (size_t i = 0; i < n; i++)
      vec[i] = GF->mod(temp[i]).lift();

    for (size_t i = n; i > 0; i--) {
      if (vec[i-1] != 0)
	{
	  R inv = GF->inverse(vec[i-1]);
	  for (size_t j = 0; j < i-1; j++)
	    vec[j] = GF->mod(GF->mul(vec[j], inv)).lift();
	  vec[i-1] = 1;
	  break;
	}
    }

    return vec;
  }

  QuadForm<T, n> get_neighbor(R t, Isometry<T, n>& s) const
    {
      Vector<R, n> vec = this->isotropic_vector(t);
        return build_neighbor(vec, s);
    }

  QuadForm<T, n> build_neighbor(Vector<R, n>& vec2, Isometry<T, n>& s) const
    {
        T p = GF->prime();
        // T pp = p*p;
        SquareMatrix<T, n> qq;

        // Convert isotropic vector into the correct domain.
        Vector<T, n> vec;
	for (size_t i = 0; i < n; i++)
	  vec[i] = GF->mod(vec2[i]).lift();

#ifdef DEBUG
        assert( q.evaluate(vec) % p == 0 );
#endif
	
	for (size_t i = 0; i < n-1; i++)
	  if (vec[i] > (p >> 1)) vec[i] -=p;

	// set s with an isometry
	// whose first column is our isotrpic vector
	size_t pivot = n-1;
	while (vec[pivot] == 0) pivot--;
#ifdef DEBUG
	assert( pivot >= 0 );
	assert( vec[pivot] == 1);
#endif
	for (size_t i = 0; i < n; i++)
	  s(i,pivot) = vec[i];
	 
	//	s.swap_cols(0, pivot);
	// need to adjust determinant for s to be in SO
	
	qq = s.transform(q.bilinear_form(), 1);

#ifdef DEBUG
        assert( qq(0,0) % p == 0 );
#endif

	// Stub
	// !! TODO - build qq to be the neighbor
	// using appropriate isometries

	QuadForm<T, n> retval(qq);
	
        if (std::is_same<T,Z64>::value)
        {
            // If we're not using arbitrary precision, throw an exception if
            // the discriminant of the p-neighbor isn't correct.
            if (retval.discriminant() != this->disc)
            {
                throw std::overflow_error(
                    "An overflow has occurred. The p-neighbor's discriminant "
                    "does not match the original.");
            }
        }
        return retval;
    }

private:
  std::shared_ptr<Fp<R,S>> GF;
  QuadForm<T, n> q;
  T disc;
  SquareMatrix< FpElement<R, S> , n> b;

  VectorFp< R, S, n> vec;

  // The 2-isotropic vectors were stored in binary within each of the
  // coordinates of `vec` and so we use this function to unpack them into
  // actual 2-isotropic vectors.
  Vector<R, n> isotropic_vector_p2(R t) const
  {
    Vector<R, n> res;
    res[n-1] = 1;
    
    // Stub
    // !! TODO - do something appropriate here
       
    return res;
  }
};

#endif // __NEIGHBOR_MANAGER_H
