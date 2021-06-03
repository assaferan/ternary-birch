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

template<typename R>
class Eigenvector;

template<typename R>
class EigenvectorContainer;

template<typename R, typename S>
class F2;

template<typename R, typename S>
class Fp;

template<typename R, typename S>
class F2Element;

template<typename R, typename S>
class FpElement;

template<typename R, size_t n>
class Genus;

template<typename R, size_t n>
class GenusRep;

template<typename Key>
class HashMap;

template<typename R, size_t n>
class Isometry;

template<typename R, typename S, typename T, size_t n>
class IsometrySequence;

template<typename R>
class Math;

template<typename R>
class Matrix;

template<typename R, typename S>
class MatrixF2;

template<typename R, typename S>
class MatrixFp;

template<typename R, typename S, typename T, size_t n>
class NeighborManager;

template<typename R, size_t n>
class QuadForm;

template<typename R, typename S, size_t n>
class QuadFormFp;

template<typename R>
class Rational;

class SetCover;

template<typename R>
class Spinor;

// This is not just a special case of Matrix<R>.
// Since the size is known, we can allocate memory statically.
// !! - TODO - can probably minimize code duplication by
// templating also the memory allocation
// alternatively, make both inherit from the same abstract class
template<typename R, size_t n>
class SquareMatrix;

template<typename R, typename S, size_t n>
class SquareMatrixFp;

template<typename R, size_t n>
class Vector;

template<typename R, typename S, size_t n>
class VectorFp;

/* Struct definitions */

template<typename R>
struct PrimeSymbol {
    R p;
    int power;
    bool ramified;
};

/* Templated class types */

// Finite fields.
typedef Fp<W16,W32>  W16_Fp;
typedef Fp<W32,W64>  W32_Fp;
typedef Fp<W64,W128> W64_Fp;
typedef F2<W16,W32>  W16_F2;

// Finite field elements
typedef FpElement<W16,W32>  W16_FpElement;
typedef FpElement<W32,W64>  W32_FpElement;
typedef FpElement<W64,W128> W64_FpElement;

// Genus
template<size_t n>
using Z_Genus = Genus<Z,n>;
template<size_t n>
using Z64_Genus = Genus<Z64,n>;
template<size_t n>
using Z128_Genus = Genus<Z128,n>;

// Genus representatives
template<size_t n>
using Z_GenusRep = GenusRep<Z, n>;
template<size_t n>
using Z64_GenusRep = GenusRep<Z64, n>;
template<size_t n>
using Z128_GenusRep = GenusRep<Z128, n>;

// Isometries.
template<size_t n>
using Z_Isometry = Isometry<Z,n>;
template<size_t n>
using Z64_Isometry = Isometry<Z64,n>;
template<size_t n>
using Z128_Isometry = Isometry<Z128,n>;

// Math.
typedef Math<Z> Z_Math;
typedef Math<Z64> Z64_Math;
typedef Math<Z128> Z128_Math;

// Variable size matrices
typedef Matrix<Z> Z_Matrix;
typedef Matrix<Z64> Z64_Matrix;
typedef Matrix<Z128> Z128_Matrix;

typedef MatrixFp< W16, W32> W16_MatrixFp;
typedef MatrixFp< W32, W64> W32_MatrixFp;
typedef MatrixFp< W64, W128> W64_MatrixFp;

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

// Prime symbols
typedef PrimeSymbol<Z>   Z_PrimeSymbol;
typedef PrimeSymbol<Z64> Z64_PrimeSymbol;
typedef PrimeSymbol<Z128> Z128_PrimeSymbol;

// Quadratic forms over the integers.
template<size_t n>
using Z_QuadForm = QuadForm<Z,n>;
template<size_t n>
using Z64_QuadForm = QuadForm<Z64,n>;
template<size_t n>
using Z128_QuadForm = QuadForm<Z128,n>;

// Quadratic forms over a finite field.
template <size_t n>
using W16_QuadForm = QuadFormFp<W16,W32,n>;
template <size_t n>
using W32_QuadForm = QuadFormFp<W32,W64,n>;
template <size_t n>
using W64_QuadForm = QuadFormFp<W64,W128,n>;

// Rational numbers
typedef Rational<Z> Z_Rational;
typedef Rational<Z64> Z64_Rational;
typedef Rational<Z128> Z128_Rational;

// Square matrices.
template<size_t n>
using W16_SquareMatrix = SquareMatrix<W16, n>;
template<size_t n>
using W32_SquareMatrix = SquareMatrix<W32, n>;
template<size_t n>
using W64_SquareMatrix = SquareMatrix<W64, n>;
template<size_t n>
using Z_SquareMatrix = SquareMatrix<Z, n>;
template<size_t n>
using Z64_SquareMatrix = SquareMatrix<Z64, n>;
template<size_t n>
using Z128_SquareMatrix = SquareMatrix<Z128, n>;

// Square matrices of Finite field elements
template<size_t n>
using W16_SquareMatrixFp = SquareMatrixFp< W16, W32, n>;
template<size_t n>
using W32_SquareMatrixFp = SquareMatrixFp< W32, W64, n>;
template<size_t n>
using W64_SquareMatrixFp = SquareMatrixFp< W64, W128, n>;

// Vectors.
template<size_t n>
using W16_Vector = Vector<W16, n>;
template<size_t n>
using W32_Vector = Vector<W32, n>;
template<size_t n>
using W64_Vector = Vector<W64, n>;
template<size_t n>
using Z_Vector = Vector<Z, n>;
template<size_t n>
using Z64_Vector = Vector<Z64, n>;
template<size_t n>
using Z128_Vector = Vector<Z128, n>;

// Vectors of Finite field elements
template<size_t n>
using W16_VectorFp = VectorFp< W16, W32, n>;
template<size_t n>
using W32_VectorFp = VectorFp< W32, W64, n>;
template<size_t n>
using W64_VectorFp = VectorFp< W64, W128, n>;

#endif // __BIRCH_H_
