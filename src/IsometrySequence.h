#ifndef __ISOMETRY_ITERATOR_H_
#define __ISOMETRY_ITERATOR_H_

#include <iterator>
#include "birch.h"

template<typename T, size_t n>
class IsometrySequenceData
{
public:
  Isometry<T,n> isometry;
    T denominator;
    size_t src;
    size_t dst;
};

template<typename R, typename S, typename T, size_t n>
class IsometrySequence
{
public:
  IsometrySequence(std::shared_ptr<Genus<T, n>> genus, const T& p)
    {
        this->prime = birch_util::convert_Integer<T,R>(p);
        this->primeT = p;

        if (this->prime == 2)
            this->GF = std::make_shared<W16_F2>(2, genus->seed());
        else
            this->GF = std::make_shared<Fp<R,S>>(prime, genus->seed(), true);

        this->genus_ = genus;

        this->current_rep = 0;
        this->current_neighbor = 0;

        const GenusRep<T, n>& cur = this->genus_->representative(this->current_rep);
        this->manager_ = std::make_shared<NeighborManager<R,S,T, n>>(cur.q, this->GF);
    }

    bool done() const
    {
        return current_rep >= this->genus_->size();
    }

  IsometrySequenceData<T, n> next()
    {
      IsometrySequenceData<T, n> isometry_data;

        // TODO: Come up with a better way to handle this case.
        if (this->done())
        {
            throw std::domain_error("No more isometries.");
        }

        // We assume that the current state is valid, and so we proceed by
        // computing the desired isometry.
        const GenusRep<T, n>& cur = this->genus_->representative(current_rep);
        GenusRep<T, n> foo = this->manager_->get_reduced_neighbor_rep(current_neighbor);
        size_t r = this->genus_->indexof(foo);
        const GenusRep<T, n>& rep = this->genus_->representative(r);

        // Set the isometry.
        isometry_data.isometry = cur.s * foo.s;
        isometry_data.isometry = isometry_data.isometry * rep.sinv;

        // Set the denominator.
        isometry_data.denominator = this->primeT;
        isometry_data.denominator *= birch_util::my_pow(cur.es);
        isometry_data.denominator *= birch_util::my_pow(rep.es);

        // Set the parent rep and the class of its neighbor.
        isometry_data.src = this->current_rep;
        isometry_data.dst = r;

        // Increment the neighbor and roll into the next neighbor, if
        // necessary.
        ++current_neighbor;
        if (current_neighbor > prime)
        {
            current_neighbor = 0;
            ++current_rep;

            if (!this->done())
            {
                // Update the neighbor manager if we've rolled over.
	      const GenusRep<T, n>& cur = this->genus_->representative(this->current_rep);
	      *this->manager_ = NeighborManager<R,S,T, n>(cur.q, this->GF);
            }
        }

        return isometry_data;
    }

private:
  std::shared_ptr<Genus<T,n>> genus_;
  std::shared_ptr<NeighborManager<R,S,T,n>> manager_;
    std::shared_ptr<Fp<R,S>> GF;
    R prime;
    T primeT;
    size_t current_rep;
    R current_neighbor;
};

#endif // __ISOMETRY_ITERATOR_H_
