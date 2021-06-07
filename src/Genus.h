#ifndef __GENUS_H_
#define __GENUS_H_

#include <unordered_map>

#include "birch.h"
#include "Eigenvector.h"
#include "HashMap.h"
#include "Isometry.h"
#include "Math.h"
#include "NeighborManager.h"
#include "Rational.h"
#include "Spinor.h"

template<typename R, size_t n>
class GenusRep
{
public:
  GenusRep() = default;
  GenusRep(const GenusRep<R, n>& genus) = default;
  GenusRep(GenusRep<R, n>&& genus) = default;

  QuadForm<R, n> q;
  Isometry<R, n> s;
  Isometry<R, n> sinv;
  Z64 parent;
  R p;
  std::map<R,int> es;
};

template<typename R, size_t n>
class Genus
{
  template<typename T, size_t N>
  friend class Genus;

public:
  // c-tors
  Genus() = default;

  Genus(const QuadForm<R, n>& q,
	const std::vector<PrimeSymbol<R>>& symbols, W64 seed=0);

  // copy c-tor
  template<typename T>
  Genus(const Genus<T, n>& src);

  template<typename T>
  static Genus<T, n> convert(const Genus<R, n>& src)
  {
    return Genus<T, n>(src);
  }

  // access
  size_t size(void) const
  {
    return this->hash->keys().size();
  }

  W64 seed(void) const
  {
    return this->seed_;
  }

  std::map<R,size_t> dimension_map(void) const
  {
    std::map<R,size_t> temp;
    size_t num_conductors = this->conductors.size();
    for (size_t k=0; k<num_conductors; k++)
      {
	temp[this->conductors[k]] = this->dims[k];
      }
    return temp;
  }

  std::map<R,std::vector<int>> hecke_matrix_dense(const R& p) const
  {
    if (this->disc % p == 0)
      {
	throw std::invalid_argument("Prime must not divide the discriminant.");
      }
    return this->hecke_matrix_dense_internal(p);
  }

  std::map<R,std::vector<std::vector<int>>> hecke_matrix_sparse(const R& p) const
  {
    if (this->disc % p == 0)
      {
	throw std::invalid_argument("Prime must not divide the discriminant.");
      }
    return this->hecke_matrix_sparse_internal(p);
  }

  Eigenvector<R> eigenvector(const std::vector<Z32>&, const R& ) const;

  std::vector<Z32> eigenvalues(EigenvectorManager<R, n>&, const R&) const;

  const GenusRep<R, n>& representative(size_t idx) const
  {
    return this->hash->get(idx);
  }

  size_t indexof(const GenusRep<R, n>& rep) const
  {
    return this->hash->indexof(rep);
  }

protected:
  R disc;
  std::vector<R> prime_divisors;
  std::vector<R> conductors;
  std::vector<size_t> dims;
  std::vector<std::vector<size_t>> num_auts;
  std::vector<std::vector<int>> lut_positions;
  Rational<Z> mass;
  std::unique_ptr<HashMap<W16>> spinor_primes;
  std::unique_ptr<HashMap<GenusRep<R, n>>> hash;
  std::unique_ptr<Spinor<R>> spinor;
  // hashing the invariants
  std::unique_ptr<HashMap<GenusRep<R, n>>> inv_hash;
  // mapping from the invariant hash to the genus representatives
  std::unordered_map< size_t, size_t > inv_map;
  W64 seed_;

  template<typename S, typename T>
  std::vector<Z32> _eigenvectors(EigenvectorManager<R, n>&,
				 std::shared_ptr<Fp<S,T>>, const R& ) const;

  // TODO: Add the actual mass formula here for reference.
  Rational<Z> get_mass(const QuadForm<R, n>&,
		       const std::vector<PrimeSymbol<R>>&);

  Rational<Z> local_factor(const Matrix< Rational<R> > & g,
			   const R & p);

  Rational<Z> combine(const QuadForm<R, n>& q,
		      const R & p);

  std::map<R,std::vector<std::vector<int>>>
  hecke_matrix_sparse_internal(const R& ) const;

  std::map<R,std::vector<int>> hecke_matrix_dense_internal(const R&) const;

  static std::set<R> witt_to_hasse(const R &,
				   const std::set<std::pair<R, int> > &);

  std::vector< Matrix<R> > decomposition(size_t k) const;
  
  std::vector< Matrix<R> >
  decomposition_recurse(const Matrix<R> & V_basis,
			const R & p, size_t k) const;
};

template<typename R, size_t n>
bool operator==(const GenusRep<R, n>& a, const GenusRep<R, n>& b)
{
    return a.q == b.q;
}

#include "Genus.inl"

#endif // __GENUS_H_
