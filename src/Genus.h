#ifndef __GENUS_H_
#define __GENUS_H_

#include "birch.h"
#include "Math.h"
#include "HashMap.h"
#include "Spinor.h"
#include "Isometry.h"
#include "NeighborManager.h"
#include "Eigenvector.h"

template<typename R, size_t Rank>
class GenusRep
{
public:
    GenusRep() = default;
    GenusRep(const GenusRep<R, Rank>& genus) = default;
    GenusRep(GenusRep<R, Rank>&& genus) = default;

    QuadForm<R, Rank> q;
    Isometry<R> s;
    Isometry<R> sinv;
    Z64 parent;
    R p;
    std::map<R,int> es;
};

template<typename R, size_t Rank>
class Genus
{
    template<typename T, size_t N>
    friend class Genus;

public:
    Genus() = default;

    Genus(const QuadForm<R, Rank>& q,
	  const std::vector<PrimeSymbol<R>>& symbols, W64 seed=0);

    template<typename T>
    Genus(const Genus<T>& src)
    {
        // Convert the discriminant.
        this->disc = birch_util::convert_Integer<T,R>(src.disc);

        // Convert the prime divisors.
        for (const T& p : src.prime_divisors)
        {
            this->prime_divisors.push_back(birch_util::convert_Integer<T,R>(p));
        }

        // Convert the conductors.
        for (const T& cond : src.conductors)
        {
            this->conductors.push_back(birch_util::convert_Integer<T,R>(cond));
        }

        // Copy dimensions.
        this->dims = src.dims;

        // Copy automorphisms counts.
        this->num_auts = src.num_auts;

        // Copy lookup table dimensions.
        this->lut_positions = src.lut_positions;

        // Copy mass.
        this->mass_x24 = src.mass_x24;

        // Build a copy of the spinor primes hash table.
        this->spinor_primes = std::unique_ptr<HashMap<W16>>(new HashMap<W16>(src.spinor_primes->size()));
        for (W16 x : src.spinor_primes->keys())
        {
            this->spinor_primes->add(x);
        }

        // Build a copy of the genus representatives hash table.
        this->hash = std::unique_ptr<HashMap<GenusRep<R>>>(new HashMap<GenusRep<R>>(src.hash->size()));
        for (const GenusRep<T>& rep : src.hash->keys())
        {
            this->hash->add(birch_util::convert_GenusRep<T,R>(rep));
        }

        // Create Spinor class.
        std::vector<R> primes;
        primes.reserve(src.spinor->primes().size());
        for (const T& p : src.spinor->primes())
        {
            primes.push_back(birch_util::convert_Integer<T,R>(p));
        }
        this->spinor = std::unique_ptr<Spinor<R>>(new Spinor<R>(primes));

        // Copy seed.
        this->seed_ = src.seed_;
    }

    template<typename T>
    static Genus<T> convert(const Genus<R>& src)
    {
        return Genus<T>(src);
    }

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

    Eigenvector<R> eigenvector(const std::vector<Z32>& vec, const R& conductor) const
    {
        size_t num_conductors = this->conductors.size();
        bool found = false;

        size_t k;
        for (k=0; k<num_conductors; k++)
        {
            if (this->conductors[k] == conductor)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            throw std::invalid_argument("Invalid conductor.");
        }

        size_t dim = this->dims[k];
        if (dim != vec.size())
        {
            throw std::invalid_argument("Eigenvector has incorrect dimension.");
        }

        size_t fulldim = this->size();

        std::vector<Z32> temp(this->size());
        const std::vector<int>& lut = this->lut_positions[k];

        for (size_t n=0; n<fulldim; n++)
        {
            if (lut[n] != -1)
            {
                temp[n] = vec[lut[n]];
            }
        }

        return Eigenvector<R>(std::move(temp), k);
    }

    std::vector<Z32> eigenvalues(EigenvectorManager<R>& vector_manager, const R& p) const
    {
        R bits16 = birch_util::convert_Integer<Z64,R>(1LL << 16);
        R bits32 = birch_util::convert_Integer<Z64,R>(1LL << 32);

        if (p == 2)
        {
            W16 prime = 2;
            std::shared_ptr<W16_F2> GF = std::make_shared<W16_F2>(prime, this->seed());
            return this->_eigenvectors<W16,W32>(vector_manager, GF, p);
        }
        else if (p < bits16)
        {
            W16 prime = birch_util::convert_Integer<R,W16>(p);
            std::shared_ptr<W16_Fp> GF = std::make_shared<W16_Fp>(prime, this->seed(), true);
            return this->_eigenvectors<W16,W32>(vector_manager, GF, p);
        }
        else if (p < bits32)
        {
            W32 prime = birch_util::convert_Integer<R,W32>(p);
            std::shared_ptr<W32_Fp> GF = std::make_shared<W32_Fp>(prime, this->seed(), false);
            return this->_eigenvectors<W32,W64>(vector_manager, GF, p);
        }
        else
        {
            W64 prime = birch_util::convert_Integer<R,W64>(p);
            std::shared_ptr<W64_Fp> GF = std::make_shared<W64_Fp>(prime, this->seed(), false);
            return this->_eigenvectors<W64,W128>(vector_manager, GF, p);
        }
    }

    const GenusRep<R>& representative(size_t n) const
    {
        return this->hash->get(n);
    }

    size_t indexof(const GenusRep<R>& rep) const
    {
        return this->hash->indexof(rep);
    }

private:
    R disc;
    std::vector<R> prime_divisors;
    std::vector<R> conductors;
    std::vector<size_t> dims;
    std::vector<std::vector<size_t>> num_auts;
    std::vector<std::vector<int>> lut_positions;
    Z mass_x24;
    std::unique_ptr<HashMap<W16>> spinor_primes;
    std::unique_ptr<HashMap<GenusRep<R>>> hash;
    std::unique_ptr<Spinor<R>> spinor;
    W64 seed_;

    template<typename S, typename T>
    std::vector<Z32> _eigenvectors(EigenvectorManager<R>& vector_manager, std::shared_ptr<Fp<S,T>> GF, const R& p) const
    {
        std::vector<Z32> eigenvalues(vector_manager.size());

        S prime = GF->prime();

        const GenusRep<R>& mother = this->hash->get(0);

        const Z32 *stride_ptr = vector_manager.strided_eigenvectors.data();

        size_t num_indices = vector_manager.indices.size();
        for (size_t index=0; index<num_indices; index++)
        {
            size_t npos = static_cast<size_t>(vector_manager.indices[index]);
            const GenusRep<R>& cur = this->hash->get(npos);
            NeighborManager<S,T,R> neighbor_manager(cur.q, GF);
            for (W64 t=0; t<=prime; t++)
            {
                GenusRep<R> foo = neighbor_manager.get_reduced_neighbor_rep((S)t);

                size_t rpos = this->hash->indexof(foo);
                size_t offset = vector_manager.stride * rpos;
                __builtin_prefetch(stride_ptr + offset, 0, 0);

                W64 spin_vals;
                if (unlikely(rpos == npos))
                {
                    spin_vals = this->spinor->norm(foo.q, foo.s, p);
                }
                else
                {
                    const GenusRep<R>& rep = this->hash->get(rpos);
                    foo.s = cur.s * foo.s;
                    R scalar = p;

                    foo.s = foo.s * rep.sinv;

                    scalar *= birch_util::my_pow(cur.es);
                    scalar *= birch_util::my_pow(rep.es);

                    spin_vals = this->spinor->norm(mother.q, foo.s, scalar);
                }

                for (Z64 vpos : vector_manager.position_lut[index])
                {
                    W64 cond = vector_manager.conductors[vpos];
                    Z32 value = birch_util::char_val(spin_vals & cond);
                    Z32 coord = vector_manager.strided_eigenvectors[offset + vpos];
                    if (likely(coord))
                    {
                        eigenvalues[vpos] += (value * coord);
                    }
                }
            }

            // Divide out the coordinate associated to the eigenvector to
            // recover the actual eigenvalue.
            for (Z64 vpos : vector_manager.position_lut[index])
            {
                size_t offset = vector_manager.stride * npos;
                Z32 coord = vector_manager.strided_eigenvectors[offset + vpos];
                assert( eigenvalues[vpos] % coord == 0 );
                eigenvalues[vpos] /= coord;
            }
        }

        return eigenvalues;
    }

    // TODO: Add the actual mass formula here for reference.
  Z get_mass(const QuadForm<R, Rank>&, const std::vector<PrimeSymbol<R>>&);

    std::map<R,std::vector<std::vector<int>>> hecke_matrix_sparse_internal(const R& p) const
    {
        size_t num_conductors = this->conductors.size();
        size_t num_primes = this->prime_divisors.size();

        std::vector<std::vector<int>> data(num_conductors);
        std::vector<std::vector<int>> indptr;
        std::vector<std::vector<int>> indices(num_conductors);

        W16 prime = birch_util::convert_Integer<R,W16>(p);

        std::shared_ptr<W16_Fp> GF;
        if (prime == 2)
            GF = std::make_shared<W16_F2>(2, this->seed());
        else
            GF = std::make_shared<W16_Fp>((W16)prime, this->seed(), true);

        std::vector<W64> all_spin_vals;
        all_spin_vals.reserve(prime+1);

        std::vector<std::vector<int>> rowdata;
        for (int dim : this->dims)
        {
            rowdata.push_back(std::vector<int>(dim));
            indptr.push_back(std::vector<int>(dim+1, 0));
        }

        const GenusRep<R>& mother = this->hash->keys()[0];
        size_t num_reps = this->size();
        for (size_t n=0; n<num_reps; n++)
        {
            const GenusRep<R>& cur = this->hash->get(n);
            NeighborManager<W16,W32,R> manager(cur.q, GF);

            for (W16 t=0; t<=prime; t++)
            {
                GenusRep<R> foo = manager.get_reduced_neighbor_rep(t);

                #ifdef DEBUG
                assert( foo.s.is_isometry(cur.q, foo.q, p*p) );
                #endif

                size_t r = this->hash->indexof(foo);

                #ifdef DEBUG
                assert( r < this->size() );
                #endif

                W64 spin_vals;
                if (r == n)
                {
                    spin_vals = this->spinor->norm(foo.q, foo.s, p);
                }
                else
                {
                    const GenusRep<R>& rep = this->hash->get(r);
                    foo.s = cur.s * foo.s;
                    R scalar = p;

                    #ifdef DEBUG
                    R temp_scalar = p*p;
                    R temp = birch_util::my_pow(cur.es);
                    temp_scalar *= temp * temp;
                    assert( foo.s.is_isometry(mother.q, foo.q, temp_scalar) );
                    #endif

                    foo.s = foo.s * rep.sinv;

                    #ifdef DEBUG
                    temp = birch_util::my_pow(rep.es);
                    temp_scalar *= temp * temp;
                    assert( foo.s.is_isometry(mother.q, mother.q, temp_scalar) );
                    #endif

                    scalar *= birch_util::my_pow(cur.es);
                    scalar *= birch_util::my_pow(rep.es);

                    #ifdef DEBUG
                    assert( scalar*scalar == temp_scalar );
                    #endif

                    spin_vals = this->spinor->norm(mother.q, foo.s, scalar);
                }

                all_spin_vals.push_back((r << num_primes) | spin_vals);
            }

            for (size_t k=0; k<num_conductors; k++)
            {
                const std::vector<int>& lut = this->lut_positions[k];
                int npos = lut[n];
                if (npos == -1) continue;

                // Populate the row data.
                std::vector<int>& row = rowdata[k];
                for (W64 x : all_spin_vals)
                {
                    int r = x >> num_primes;
                    int rpos = lut[r];
                    if (rpos == -1) continue;

                    int value = birch_util::char_val(x & k);
                    row[rpos] += value;
                }

                // Update data and indices with the nonzero values.
                size_t nnz = 0;
                size_t pos = 0;
                std::vector<int>& data_k = data[k];
                std::vector<int>& indices_k = indices[k];
                for (int x : row)
                {
                    if (x)
                    {
                        data_k.push_back(x);
                        indices_k.push_back(pos);
                        row[pos] = 0; // Clear the nonzero entry.
                        ++nnz;
                    }
                    ++pos;
                }

                // Update indptr
                indptr[k][npos+1] = indptr[k][npos] + nnz;
            }

            all_spin_vals.clear();
        }

        std::map<R,std::vector<std::vector<int>>> csr_matrices;
        for (size_t k=0; k<num_conductors; k++)
        {
            const R& cond = this->conductors[k];
            csr_matrices[cond] = std::vector<std::vector<int>>();
            csr_matrices[cond].push_back(data[k]);
            csr_matrices[cond].push_back(indices[k]);
            csr_matrices[cond].push_back(indptr[k]);
        }
        return csr_matrices;
    }

    std::map<R,std::vector<int>> hecke_matrix_dense_internal(const R& p) const
    {
        size_t num_conductors = this->conductors.size();
        size_t num_primes = this->prime_divisors.size();

        // Allocate memory for the Hecke matrices and create a vector to store
        // pointers to the raw matrix data.
        std::vector<int*> hecke_ptr;
        hecke_ptr.reserve(num_conductors);
        std::vector<std::vector<int>> hecke_matrices;
        for (size_t k=0; k<num_conductors; k++)
        {
            size_t dim = this->dims[k];
            hecke_matrices.push_back(std::vector<int>(dim * dim));
            hecke_ptr.push_back(hecke_matrices.back().data());
        }

        W16 prime = birch_util::convert_Integer<R,W16>(p);
        std::vector<W64> all_spin_vals;
        all_spin_vals.reserve(prime+1);

        std::shared_ptr<W16_Fp> GF;
        if (prime == 2)
            GF = std::make_shared<W16_F2>(2, this->seed());
        else
            GF = std::make_shared<W16_Fp>((W16)prime, this->seed(), true);

        const GenusRep<R>& mother = this->hash->keys()[0];
        size_t num_reps = this->size();

        // Create hash tables for storing isotropic vectors to be skipped
        // at later iterations.
        std::vector<HashMap<W16_Vector3>> vector_hash(num_reps);

        for (size_t n=0; n<num_reps; n++)
        {
            const GenusRep<R>& cur = this->hash->get(n);
            NeighborManager<W16,W32,R> manager(cur.q, GF);

            for (W16 t=0; t<=prime; t++)
            {
                GenusRep<R> foo;
                W16_Vector3 vec = manager.isotropic_vector(t);
                vec.x = GF->mod(vec.x);
                vec.y = GF->mod(vec.y);
                vec.z = GF->mod(vec.z);

                // If this vector has already been identified, skip it. This
                // avoids unnecessarily computing the neighbor, reducing it,
                // and testing for isometry. The Hermitian symmetry property
                // of the Hecke matrix will account for this once we finish
                // processing neighbors.
                if (vector_hash[n].exists(vec)) continue;

                // Build the neighbor and reduce it.
                foo.q = manager.build_neighbor(vec, foo.s);
                foo.q = QuadForm<R>::reduce(foo.q, foo.s);

                #ifdef DEBUG
                assert( foo.s.is_isometry(cur.q, foo.q, p*p) );
                #endif

                size_t r = this->hash->indexof(foo);

                #ifdef DEBUG
                assert( r < this->size() );
                #endif

                W64 spin_vals;
                if (r > n)
                {
                    W16_Vector3 result = manager.transform_vector(foo, vec);
                    vector_hash[r].add(result);

                    const GenusRep<R>& rep = this->hash->get(r);
                    foo.s = cur.s * foo.s;
                    R scalar = p;

                    #ifdef DEBUG
                    R temp_scalar = p*p;
                    R temp = birch_util::my_pow(cur.es);
                    temp_scalar *= temp * temp;
                    assert( foo.s.is_isometry(mother.q, foo.q, temp_scalar) );
                    #endif

                    foo.s = foo.s * rep.sinv;

                    #ifdef DEBUG
                    temp = birch_util::my_pow(rep.es);
                    temp_scalar *= temp * temp;
                    assert( foo.s.is_isometry(mother.q, mother.q, temp_scalar) );
                    #endif

                    scalar *= birch_util::my_pow(cur.es);
                    scalar *= birch_util::my_pow(rep.es);

                    #ifdef DEBUG
                    assert( scalar*scalar == temp_scalar );
                    #endif

                    spin_vals = this->spinor->norm(mother.q, foo.s, scalar);
                }
                else if (r == n)
                {
                    spin_vals = this->spinor->norm(foo.q, foo.s, p);
                }
                else continue;

                all_spin_vals.push_back((r << num_primes) | spin_vals);
            }

            for (size_t k=0; k<num_conductors; k++)
            {
                const std::vector<int>& lut = this->lut_positions[k];
                int npos = lut[n];
                if (unlikely(npos == -1)) continue;

                int *row = hecke_ptr[k];

                for (W64 x : all_spin_vals)
                {
                    int r = x >> num_primes;
                    int rpos = lut[r];
                    if (unlikely(rpos == -1)) continue;

                    row[rpos] += birch_util::char_val(x & k);
                }

                hecke_ptr[k] += this->dims[k];
            }

            all_spin_vals.clear();
        }

        // Copy the upper diagonal entries to the lower diagonal using the
        // Hermitian symmetry property and then move the matrix into an
        // associatively map before returning.
        std::map<R,std::vector<int>> matrices;
        for (size_t k=0; k<num_conductors; k++)
        {
            std::vector<int>& matrix = hecke_matrices[k];
            size_t dim = this->dims[k];
            size_t dim2 = dim * dim;
            const std::vector<size_t>& auts = this->num_auts[k];

            // Copy upper diagonal matrix to the lower diagonal.
            for (size_t start=0, row=0; start<dim2; start+=dim+1, row++)
            {
                int row_auts = auts[row];
                for (size_t dst=start+dim, src=start+1, col=row+1; col<dim; src++, col++, dst+=dim)
                {
                    if (matrix[src])
                    {
                        int col_auts = auts[col];
                        if (col_auts == row_auts)
                        {
                            matrix[dst] = matrix[src];
                        }
                        else
                        {
                            #ifdef DEBUG
                            assert( (matrix[src] * col_auts) % row_auts == 0 );
                            #endif

                            matrix[dst] = matrix[src] * col_auts / row_auts;
                        }
                    }
                }
            }

            // Move the matrix in the corresponding entry in the map.
            matrices[this->conductors[k]] = std::move(hecke_matrices[k]);
        }
        return matrices;
    }
};

template<typename R>
bool operator==(const GenusRep<R>& a, const GenusRep<R>& b)
{
    return a.q == b.q;
}

#include "Genus.inl"

#endif // __GENUS_H_
