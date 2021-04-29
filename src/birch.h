#ifndef __BIRCH_H_
#define __BIRCH_H_

#include <stdint.h>
#include <iostream>
#include <chrono>
#include <cassert>
#include <vector>
#include <map>
#include <array>
#include <random>
#include <memory>
#include <gmpxx.h>

/* Type definitions */

typedef mpz_class W;
typedef uint16_t W16;
typedef uint32_t W32;
typedef uint64_t W64;
typedef __uint128_t W128;

typedef mpz_class Z;
typedef int16_t Z16;
typedef int32_t Z32;
typedef int64_t Z64;
typedef __int128_t Z128;

/* Builtins */

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

// Constants
constexpr W64 FNV_OFFSET = 0x811c9dc5;
constexpr W64 FNV_PRIME  = 0x01000193;

/* Class declarations */

class SetCover;

template<typename R>
class Math;

template<typename R, size_t n>
class Isometry;

template<typename R, size_t n>
class QuadForm;

template<typename R, typename S, size_t n>
class QuadFormFp;

template<typename R, typename S>
class Fp;

template<typename R, typename S>
class F2;

template<typename R>
class Eigenvector;

template<typename R>
class EigenvectorContainer;

template<typename R, typename S, typename T, size_t n>
class NeighborManager;

template<typename Key>
class HashMap;

template<typename R, size_t n>
class Genus;

template<typename R>
class Spinor;

template<typename R, size_t n>
class GenusRep;

template<typename R, typename S, typename T, size_t n>
class IsometrySequence;



/* Struct definitions */

template<typename R>
struct Vector3 {
    R x; R y; R z;
};

template<typename R, size_t n>
struct Vector {
  R v[n];
  // backward compatibility
  R x; R y; R z;
  const R& operator[](size_t i) const {return v[i]; }
  R& operator[](size_t i) {return v[i];}
  bool operator==(const Vector<R,n> & other) const
  {
    for (size_t i = 0; i < n; i++)
      if (this->v[i] != other[i]) return false;
    return true;
  }
};

namespace std
{
  template<typename R, size_t n>
  struct hash<Vector<R, n> >
    {
      Z64 operator()(const Vector<R,n>& vec) const
        {
            Z64 fnv = FNV_OFFSET;
	    for (size_t i = 0; i < n; i++)
	      fnv = (fnv ^ vec.v[i]) * FNV_PRIME;
            
            return fnv;
        }
    };
}

template<typename R>
struct PrimeSymbol {
    R p;
    int power;
    bool ramified;
};

/* Templated class types */

// Isometries.
template<size_t n>
using Z_Isometry = Isometry<Z,n>;

template<size_t n>
using Z64_Isometry = Isometry<Z64,n>;

// Quadratic forms over the integers.
template<size_t n>
using Z_QuadForm = QuadForm<Z,n>;

template<size_t n>
using Z64_QuadForm = QuadForm<Z,n>;

// Quadratic forms over a finite field.
template <size_t n>
using W16_QuadForm = QuadFormFp<W16,W32,n>;
template <size_t n>
using W32_QuadForm = QuadFormFp<W32,W64,n>;
template <size_t n>
using W64_QuadForm = QuadFormFp<W64,W128,n>;

// Vectors.
template<size_t n>
using W16_Vector = Vector<W16, n>;
template<size_t n>
using W32_Vector = Vector<W32, n>;
template<size_t n>
using W64_Vector = Vector<W64, n>;
template<size_t n>
using Z_Vector = Vector<Z, n>;

// Finite fields.
typedef Fp<W16,W32>  W16_Fp;
typedef Fp<W32,W64>  W32_Fp;
typedef Fp<W64,W128> W64_Fp;
typedef F2<W16,W32>  W16_F2;

// Prime symbols
typedef PrimeSymbol<Z>   Z_PrimeSymbol;
typedef PrimeSymbol<Z64> Z64_PrimeSymbol;

// Math.
typedef Math<Z> Z_Math;
typedef Math<Z64> Z64_Math;

// Neighbor managers.
template <size_t n>
using Z_W16_NeighborManager = NeighborManager<W16,W32,Z,n>;
template <size_t n>
using Z_W32_NeighborManager = NeighborManager<W32,W64,Z,n>;
template <size_t n>
using Z_W64_NeighborManager =  NeighborManager<W64,W128,Z,n>;
template <size_t n>
using Z64_W16_NeighborManager = NeighborManager<W16,W32,Z64,n>;
template <size_t n>
using Z64_W32_NeighborManager = NeighborManager<W32,W64,Z64,n>;
template <size_t n>
using Z64_W64_NeighborManager = NeighborManager<W64,W128,Z64,n>;

// Genus
template<size_t n>
using Z_Genus = Genus<Z,n>;

template<size_t n>
using Z64_Genus = Genus<Z64,n>;

// Genus representatives
template<size_t n>
using Z_GenusRep = GenusRep<Z, n>;

template<size_t n>
using Z64_GenusRep = GenusRep<Z64, n>;


#endif // __BIRCH_H_
