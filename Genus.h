#ifndef __GENUS_H_
#define __GENUS_H_

#include <memory>
#include <iomanip>
#include <cassert>
#include <unordered_set>
#include <set>
#include <map>
#include <queue>
#include <thread>
#include <mutex>
#include <chrono>
#include "QuadForm.h"
#include "NeighborIterator.h"
#include "Character.h"

template<typename R, typename F>
class GenusRep
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;

    GenusRep(QuadFormPtr qPtr);

    int64_t pos(void) const;
    void pos(int64_t x);
    QuadFormPtr quad_form(void) const;

    void set_dimension(const Character<R,F>& chi, const mpz_class& dim);
    mpz_class dimension(const Character<R,F>& chi) const;

    bool operator==(const GenusRep<R,F>& rep) const;
    bool operator<(const GenusRep<R,F>& rep) const;
private:
    QuadFormPtr q_;
    int64_t pos_ = -1;
    std::map<Character<R,F>, mpz_class> dimensionMap_;
};

template<typename R, typename F>
class Genus
{
public:
    typedef std::shared_ptr<QuadForm<R,F>> QuadFormPtr;

    Genus(const QuadForm<R,F>& q, int64_t numThreads=0);

    QuadFormPtr quad_form(void) const;

    R smallest_good_prime(void) const;

    void add_character(const Character<R,F>& chi);
    void add_genus_rep(QuadFormPtr q, QuadFormPtr neighbor=nullptr);

    void set_dimensions(GenusRep<R,F>& rep);

    void print(void) const;

    size_t size(void) const;

    mpz_class dimension(const Character<R,F>& chi) const;

    /* Computes the full genus. */
    // TODO: MAKE THIS PRIVATE. THIS IS PUBLIC FOR TESTING PURPSOES ONLY.
    void compute_genus(void);

private:
    /* The number of threads to use. */
    int64_t numThreads_;

    /* This function is passed to threads, pulls a genus representative off of
     * a queue, then computes its p-neighbors. */
    void threaded_compute_genus(const R& p, int64_t threadid);

    /* A vector of booleans used to determine whether a thread is blocked due
     * to there being no objects in the genusQueue_. */
    std::vector<bool> threadBlocked_;

    /* An unordered set which stores distinct genus representatives. */
    std::unordered_set<GenusRep<R,F>> genusReps_;

    /* A vector of shared QuadForm pointers. An unordered list of the genus. */
    std::vector<QuadFormPtr> genusVec_;

    /* A queue of shared QuadForm pointers used for multithreading. */
    std::queue<QuadFormPtr> genusQueue_;

    /* A mutex used for multithreading. */
    std::mutex genusMutex_;

    /* A pointer to the originating form. */
    QuadFormPtr q_;

    /* The discriminant. */
    R disc_;

    /* Flag indicating whether the genus has been computed in full. */
    bool computed_ = false;

    /* A set of characters to compute. */
    std::set<Character<R,F>> charSet_;

    /* A set of prime characters to compute. When representation
     * values are compuated for characters in charSet_, these will first be
     * computed and then multiplied together to obtain the desired
     * representation value. */
    std::set<Character<R,F>> primeCharSet_;
};

template<typename R, typename F>
GenusRep<R,F>::GenusRep(QuadFormPtr q)
{
    this->q_ = q;
}

template<typename R, typename F>
int64_t GenusRep<R,F>::pos(void) const
{
    return this->pos_;
}

template<typename R, typename F>
void GenusRep<R,F>::pos(int64_t x)
{
    this->pos_ = x;
}

template<typename R, typename F>
std::shared_ptr<QuadForm<R,F>> GenusRep<R,F>::quad_form(void) const
{
    return this->q_;
}

template<typename R, typename F>
bool GenusRep<R,F>::operator==(const GenusRep<R,F>& rep) const
{
    return *(rep.quad_form()) == *this->q_;
}

template<typename R, typename F>
bool GenusRep<R,F>::operator<(const GenusRep<R,F>& rep) const
{
    return *this->q_ < *(rep.quad_form());
}

template<typename R, typename F>
mpz_class GenusRep<R,F>::dimension(const Character<R,F>& chi) const
{
    if (this->dimensionMap_.count(chi) > 0)
    {
        return this->dimensionMap_.find(chi)->second;
    }
    else
    {
        return 0;
    }
}

namespace std
{
    template<typename R, typename F>
    struct hash<GenusRep<R,F>>
    {
        int64_t operator()(const GenusRep<R,F>& rep) const
        {
            return rep.quad_form()->hash_value();
        }
    };
}

template<typename R, typename F>
Genus<R,F>::Genus(const QuadForm<R,F>& q, int64_t numThreads)
{
    this->q_ = QuadForm<R,F>::reduce(q, false);
    this->disc_ = this->q_->discriminant();
    this->numThreads_ = numThreads;
}

template<typename R, typename F>
std::shared_ptr<QuadForm<R,F>> Genus<R,F>::quad_form(void) const
{
    return this->q_;
}

template<typename R, typename F>
size_t Genus<R,F>::size(void) const
{
    return this->genusSet_.size();
}

template<typename R, typename F>
void Genus<R,F>::threaded_compute_genus(const R& p, int64_t threadid)
{
    do
    {
        // Obtain the mutex.
        this->genusMutex_.lock();

        // Once the lock is obtained, indicate that the thread is not blocked.
        // This will ensure that other threads will not terminate if they
        // obtain the lock.
        this->threadBlocked_[threadid] = false;

        // Check the size of the queue to determine whether we are blocked.
        if (this->genusQueue_.size() == 0)
        {
            // If there is nothing in the queue, we're blocked.
            this->threadBlocked_[threadid] = true;

            // Determine whether all threads are blocked.
            bool done = true;
            for (bool value : this->threadBlocked_)
            {
                done = done && value;
            }

            // If all threads are blocked, release the lock and terminate.
            if (done) 
            {
                this->genusMutex_.unlock();
                return;
            }

            // Release the lock and sleep for a random number of milliseconds.
            this->genusMutex_.unlock();
            std::default_random_engine dre(threadid);
            std::uniform_int_distribution<int64_t> id(0, 20);
            std::this_thread::sleep_for(std::chrono::milliseconds(id(dre)));
            continue;
        }
        else
        {
            // Pop the next genus representative off the queue, and unlock the
            // mutex.
            QuadFormPtr cur = this->genusQueue_.front();
            this->genusQueue_.pop();
            this->genusMutex_.unlock();

            // The neighbor iterator.
            NeighborIterator<R,F> it(cur, p);

            // The first neighbor.
            QuadFormPtr pn = it.next_neighbor();
#ifdef DEBUG
            assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif

            // Loop over all p-neighbors.
            while (pn != nullptr)
            {
#ifdef DEBUG
                assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif

                // Reduce the p-neighbor.
                QuadFormPtr qq = QuadForm<R,F>::reduce(*pn);

                // Obtain the lock, add this genus representative to the queue
                // and release the lock.
                this->genusMutex_.lock();
                this->add_genus_rep(qq, pn);
                this->genusMutex_.unlock();

                // Get the next p-neighbor.
                pn = it.next_neighbor();
            }
        }
    }
    while (true);
}

template<typename R, typename F>
void Genus<R,F>::compute_genus(void)
{
    // Do nothing if the genus has already been computed.
    if (this->computed_) { return; }

    // Add the defining quadratic form as the initial genus rep.
    this->add_genus_rep(this->q_);

    // Determine the smallest good prime.
    R p = this->smallest_good_prime();

    // Should we use threads?
    if (this->numThreads_ > 0)
    {
        // Initialize a vector of booleans initialized to false. Each thread
        // will modify this vector to indicate that it is blocked on the
        // queue. Once all threads have indicated that they are blocked, each
        // thread will terminate execution.
        this->threadBlocked_ = std::vector<bool>(this->numThreads_, false);

        // Build the threads.
        std::vector<std::thread> threads;
        for (int64_t n = 0; n < this->numThreads_; n++)
        {
            threads.push_back(std::thread(&Genus<R,F>::threaded_compute_genus, this, p, n));
            this->threadBlocked_[n] = false;
        }

        // Wait for all threads to terminate.
        for (auto& t : threads)
        {
            t.join();
        }
    }
    else
    {
        // Loop over all genus representatives as the genus is built.
        int64_t index = 0;
        while (index < this->genusVec_.size())
        {
            // The current genus representative.
            QuadFormPtr cur = this->genusVec_[index++];
    
            // The neighbor iterator.
            NeighborIterator<R,F> it(cur, p);
    
            // The first neighbor.
            QuadFormPtr pn = it.next_neighbor();
    
#ifdef DEBUG
            assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif
    
            // Loop over all p-neighbors.
            while (pn != nullptr)
            {
#ifdef DEBUG
                assert( pn->isometry()->is_isometry(*this->q_, *pn) );
#endif
    
                // Reduce the p-neighbor.
                QuadFormPtr qq = QuadForm<R,F>::reduce(*pn);
    
                // Add this reduce form to the genus, if necessary.
                this->add_genus_rep(qq, pn);
    
                // Get the next p-neighbor.
                pn = it.next_neighbor();
            }
        }
    }
}

template<typename R, typename F>
void Genus<R,F>::add_genus_rep(QuadFormPtr q, QuadFormPtr neighbor)
{
    // Create a GenusRep object from the shared pointer.
    GenusRep<R,F> rep(q);

    // Is this quadratic form already in the genus? If not, add it.
    if (this->genusReps_.count(rep) == 0)
    {
        // Compose to obtain a global isometry between this genus rep and
        // the original quadratic form.
        if (neighbor != nullptr)
        {
            q->isometry()->multiply_on_left_by(neighbor->isometry());

#ifdef DEBUG
            assert( q->isometry()->is_isometry(*this->q_, *q) );
#endif
        }

        // Set the dimensions for this genus representative.
        this->set_dimensions(rep);

        // Update data structures.
        this->genusReps_.insert(rep);
        this->genusVec_.push_back(q);
        this->genusQueue_.push(q);
    }
}

template<typename R, typename F>
void GenusRep<R,F>::set_dimension(const Character<R,F>& chi, const mpz_class& dim)
{
    this->dimensionMap_[chi] = dim;
}

template<typename R, typename F>
void Genus<R,F>::add_character(const Character<R,F>& chi)
{
    // If this is a new character, add it to the set.
    if (this->charSet_.count(chi) == 0)
    {
        this->charSet_.insert(chi);

        // The list of primes associated to this character.
        std::vector<R> ps = chi.primes();

        // Loop over all of these primes and add any new primes to the prime
        // set.
        for (auto& p : ps)
        {
            Character<R,F> temp(p);
            if (this->primeCharSet_.count(temp) == 0)
            {
                this->primeCharSet_.insert(temp);
            }
        }
    }
}

template<typename R, typename F>
mpz_class Genus<R,F>::dimension(const Character<R,F>& chi) const
{
    mpz_class dim = 0;
    for (auto& rep : this->genusReps_)
    {
        dim += rep.dimension(chi);
    }
    return dim;
}

template<typename R, typename F>
void Genus<R,F>::print(void) const
{
    for (auto& chi : this->charSet_)
    {
        std::cout << std::setw(10) << chi.conductor() << "   ";
        std::cout << this->dimension(chi) << std::endl;
    }
    //for (auto& q : this->genusVec_)
    //{
    //    std::cout << q << std::endl;
    //}
}

#endif // __GENUS_H_
