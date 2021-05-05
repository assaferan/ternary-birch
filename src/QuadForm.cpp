#include <sstream>
#include <type_traits>

#include "birch.h"
#include "Isometry.h"
#include "QuadForm.h"
#include "Math.h"

template class QuadForm_Base<Z, 5>;
template class QuadFor_Base<Z64, 5>;

template<typename R, size_t n>
std::vector< QuadForm<R, 5> >
QuadForm_Base<R, n>::nipp_to_forms(NippEntry entry)
{
  std::vector< QuadForm<R, 5> > forms;
  size_t triangular[5];
  for (size_t j = 0; j < 5; j++)
    triangular[j] = j*(j-1)/2;
  typename QuadForm<R,5>::SymVec form;
  for (LatticeRecord lat : entry.lattices)
    {
      size_t form_idx = 0;
      for (size_t col = 0; col < 5; col++)
	{
	  for (size_t row = 0; row < col; row++)
	    {
	      form[form_idx++] = lat.form[5+triangular[col]+row]; 
	    }
	  form[form_idx++] = 2*lat.form[col];
	}
      forms.push_back(form);
    }
  return forms;
}

template<typename R, size_t n>
std::vector<std::vector< QuadForm<R, 5> > >
QuadForm_Base<R,n>::get_quinary_forms(const R & disc)
{
  std::vector< std::vector< QuadForm<R, 5> > > all_forms;

  std::vector<R> nipp_maxs = {0,256,270,300,322,345,400,440,480,500,513};
  size_t table_idx = 0;
  while (nipp_maxs[table_idx+1] < disc) table_idx++;
  std::ostringstream nipp_fname;
  nipp_fname << "lattice_db/nipp" << nipp_maxs[table_idx]+1 << "-";
  nipp_fname << nipp_maxs[table_idx+1] << ".txt";
  std::cout << "nipp_fname = " << nipp_fname.str() << std::endl;

  std::vector<NippEntry> nipps =
    ParseNipp::parseDisc(nipp_fname.str(), disc);
  
  for (NippEntry nipp : nipps)
    {
#ifdef DEBUG
      std::cerr << "disc = " << nipp.disc << std::endl;
      std::cerr << "genus = " << nipp.genus << std::endl;
      std::cerr << "mass = " << nipp.mass[0] << "/" << nipp.mass[1] << std::endl;
      std::cerr << "Hasse symbols = ";
      for (short int symb : nipp.HasseSymb)
	std::cerr << symb << " ";
      std::cerr << std::endl;
      std::cerr << "lattices = " << std::endl;
      for (LatticeRecord lat : nipp.lattices)
	{
	  for (Z num : lat.form)
	    std::cerr << num << " ";
	  std::cerr << ";\t" << lat.numAut << std::endl; 
	}
#endif
      all_forms.push_back(QuadForm<R, 5>::nipp_to_forms(nipp));
    }
  
  return all_forms;
}

static W64 sign_vector(const Z& x, const Z& det,
		       const std::vector<Z_PrimeSymbol>& primes)
{
    W64 vec = (x == -1);
    for (const Z_PrimeSymbol& symb : primes)
    {
        int value = Z_Math::hilbert_symbol(x, -det, symb.p);
        vec = (vec << 1) | (value == -1);
    }
    return vec;
}

// A naive GF(2) solver that looks for solutions to a specific linear equation
// from a specified starting point.
static W64 GF2_solve_naive(const std::vector<W64>& vecs, W64 start, W64 target)
{
    W64 upper = 1LL << vecs.size();
    size_t num_vecs = vecs.size();
    for (W64 i=start+1; i<upper; i++)
    {
        W64 x = 0;
        W64 mask = upper >> 1;
        for (size_t j=0; j<num_vecs; j++)
        {
            if (i & mask) x ^= vecs[j];
            mask >>= 1;
        }

        if (x == target) return i;
    }

    return 0;
}

// A simple brute force search for p^2-isotropic vectors. This can probably be
// rewritten to avoid an exhaustive search, but since we expect the primes to
// be small, this should work for now.
template<size_t n>
static Z_Vector<n> Z_isotropic_mod_pp(const Z_QuadForm<n>& q, const Z& p)
{
    Z pp = p*p;
    Z_Vector<n> vec;
    vec[n-1] = 1;

    // Try (0 0 ... 0 1) first.
    if (q.evaluate(vec) % pp == 0)
        return vec;

    // Try (0 0 .... 0 1 x) next.
    vec[n-2] = 1;
    for (Z x = 0; x < pp; x++) {
      vec[n-1] = x;
      if (q.evaluate(vec) % pp == 0)
	return vec;
    }
    
    // Lastly, try (0 0 .... 0 1 x y).
    vec[n-3] = 1;
    for (Z x = 0; x < pp; x++) {
      vec[n-2] = x;
      for (Z y = 0; y < pp; y++) {
	vec[n-1] = y;
	if (q.evaluate(vec) % pp == 0)
	  return vec;
      }
    }
    // Otherwise, return the zero vector to indicate no solution found.
    for (size_t i = 0; i < n; i++) vec[i] = 0;
    return vec;
}

template<>
Z_QuadForm<3> Z_QuadForm<3>::get_quad_form(const std::vector<Z_PrimeSymbol>& input)
{
    Z det = 1;
    Z disc = 1;
    bool two_specified = false;

    int num_ramified = 0;
    for (const Z_PrimeSymbol& symb : input)
    {
        if (symb.power % 2 == 0)
        {
            throw std::invalid_argument("Prime powers must be odd.");
        }

        if (symb.p == 2)
        {
            two_specified = true;
        }

        if (symb.ramified)
        {
            ++num_ramified;
        }

        det *= symb.p;
        for (int k=0; k<symb.power; k++)
        {
            disc *= symb.p;
        }
    }

    if (num_ramified % 2 == 0)
    {
        throw std::invalid_argument("Must specify an odd number of ramified primes.");
    }

    // Make a copy of the prime symbols to include 2 if it wasn't already
    // included. We also build the target vector encoding the desired
    // ramification properties.
    std::vector<Z_PrimeSymbol> primes;
    W64 target = 1; // Ramified at the infinite place.

    if (!two_specified)
    {
        Z_PrimeSymbol s;
        s.p = 2;
        s.power = 0;
        s.ramified = 0;
        primes.push_back(s);
        target <<= 1;
    }
    for (const Z_PrimeSymbol& symb : input)
    {
        primes.push_back(symb);
        target <<= 1;
        if (symb.ramified)
        {
            target |= 1;
        }
    }

    // The list of primes used in our factor baes.
    std::vector<Z> fullbase;

    // Create an std::vector consisting of GF(2)-vectors encoding Hilbert
    // symbol values.
    std::vector<W64> signs;

    // Add the relation for the infinite prime.
    fullbase.push_back(-1);
    signs.push_back(sign_vector(-1, det, primes));

    for (const Z_PrimeSymbol& symb : primes)
    {
        signs.push_back(sign_vector(symb.p, det, primes));
        fullbase.push_back(symb.p);
    }

    Z p = 2;
    W64 solution;
    bool done = false;
    bool added_to_end = false;

    while (!done)
    {
        solution = 0;
        do
        {
            solution = GF2_solve_naive(signs, solution, target);

            if (solution)
            {
                W64 mask = 1LL << fullbase.size();
                Z b = -1;
                Z det2 = det;

                // Construct the quadratic space associated to the solution
                // we've found.
                for (const Z& q : fullbase)
                {
                    mask >>= 1;
                    if (solution & mask)
                    {
                        b *= q;
                        if (q > 0)
                        {
                            if (det2 % q == 0) det2 /= q;
                            else
                            {
                                det2 *= q;
                            }
                        }
                    }
                }

                // Verify that the Hilbert symbols are correct at inf, 2, and
                // the primes dividing the discriminant.
                mask = 1LL << primes.size();
                for (const Z_PrimeSymbol& symb : primes)
                {
                    mask >>= 1;
                    int sign = (target & mask) ? -1 : 1;
                    if (Z_Math::hilbert_symbol(-b, -disc, symb.p) != sign)
                    {
                        throw std::runtime_error("Hilbert symbols do not check out. How did this happen?");
                    }
                }

                // We assume that the quadratic space is good outside of the
                // initial specifications. If we detect that one of our
                // auxilliary primes is ramified, we immediately flag this
                // space as not good and try another solution.
                int good = true;
                for (size_t n=primes.size()+1; n<fullbase.size(); n++)
                {
                    int sign = Z_Math::hilbert_symbol(-b, -disc, fullbase[n]);
                    if (sign == -1) good = false;
                }

                // If good, we specify as such and break out of the while loop
                // searching for solutions.
                if (good)
                {
                    done = true;
                    break;
                }
            }
        }
        while (solution != 0);

        // If we flagged done within our while-loop above, i.e. we found a
        // quadratic space with the correct ramification behavior, then we
        // break out of the factor base building loop.
        if (done) break;

        // Get the next prime not dividing the discriminant...
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        while (disc % p == 0)
        {
            mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        }

        // If we've added a prime to the end of our factor base list, we
        // remove it if it didn't lead to a solution. This helps avoid the
        // factor base growing too large, although it can lead to the size of
        // the primes in the factor base becoming larger than we'd like. (This
        // does not appear to be a problem in practice.)
        if (added_to_end)
        {
            signs.pop_back();
            fullbase.pop_back();
        }

        // ...and push it's sign vector onto the list.
        signs.push_back(sign_vector(p, det, primes));
        fullbase.push_back(p);
        added_to_end = true;
    }

    // Construct the coefficients of the diagonalized quadratic space with
    // the desired signs.
    Z b = -1;
    W64 mask = 1LL << fullbase.size();

    // Set up the base primes. These are primes which do not divide the
    // discriminant, but were included (as squares) in the discriminant of
    // the resulting quadratic space. We also include 2 in this list, as
    // additional powers of 2 can sometimes occur, even if 2 divides the
    // discriminant.
    std::vector<Z> base;
    base.push_back(2);

    for (const Z& p : fullbase)
    {
        mask >>= 1;
        if (solution & mask)
        {
            b *= p;
            if (p > 0)
            {
                if (det % p == 0) det /= p;
                else
                {
                    // Add p to the base primes.
                    det *= p;
                    base.push_back(p);
                }
            }
        }
    }

    // Verify that the Hilbert symbols are correct.
    mask = 1LL << primes.size();
    for (const Z_PrimeSymbol& symb : primes)
    {
        mask >>= 1;
        int sign = (target & mask) ? -1 : 1;
        if (Z_Math::hilbert_symbol(-b, -disc, symb.p) != sign)
        {
            throw std::runtime_error("Hilbert symbols do not check out. How did this happen?");
        }
    }

    // Start by building a diagonal form.
    Z_QuadForm<3>::SymVec form;
    Z a = 1;
    Z c = det;
    Z f = 0;
    Z g = 0;
    Z h = 0;
    form[0] = 2*a;
    form[1] = h;
    form[2] = 2*b;
    form[3] = f;
    form[4] = g;
    form[5] = 2*c;
    Z_QuadForm<3> q(form);
    Z N = q.discriminant();

    // Remove the prime squares from the discriminant for those primes not
    // dividing the intended discriminant.
    for (const Z& p : base)
    {
        Z pp = p*p;

        while (N % pp == 0)
        {
            // Find an isotropic vector mod p^2, make it a basis vector,
            // and then divide p^2 out of the discriminant.
	  Z_Vector<3> vec = Z_isotropic_mod_pp(q, p);

            #ifdef DEBUG
            assert( q.evaluate(vec) % pp == 0 );
            #endif

            if (vec[0] == 0 && vec[1] == 0 && vec[2] == 0) break;
            else if (vec[0] == 0 && vec[1] == 0)
            {
                #ifdef DEBUG
                assert( vec[2] == 1 );
                assert( g % p == 0 );
                assert( f % p == 0 );
                assert( c % pp == 0 );
                #endif

                c /= pp;
                f /= p;
                g /= p;
            }
            else if (vec[0] == 0)
            {
                b += (c*vec[2]*vec[2] + f*vec[2]);
                f += (2*c*vec[2]);
                h += (g*vec[2]);

                #ifdef DEBUG
                assert( vec[1] == 1 );
                assert( b % pp == 0 );
                assert( f % p == 0 );
                assert( h % p == 0 );
                #endif

                b /= pp;
                f /= p;
                h /= p;
            }
            else
            {
                a += (b*vec[1]*vec[1] + c*vec[2]*vec[2] + f*vec[1]*vec[2] + g*vec[2] + h*vec[1]);
                g += (2*c*vec[2] + f*vec[1]);
                h += (2*b*vec[1] + f*vec[2]);

                #ifdef DEBUG
                assert( vec[0] == 1 );
                assert( a % pp == 0 );
                assert( g % p == 0 );
                assert( h % p == 0 );
                #endif

                a /= pp;
                g /= p;
                h /= p;
            }
	    form[0] = 2*a;
	    form[1] = h;
	    form[2] = 2*b;
	    form[3] = f;
	    form[4] = g;
	    form[5] = 2*c;
            q = Z_QuadForm<3>(form);
            N = q.discriminant();
        }
    }

    // Reduce the resulting quadratic form, scale up to the correct
    // discriminant, then reduce again. The resulting form will have the
    // correct local behavior as well as the correct discriminant.
    Z_Isometry<3> s;
    q = Z_QuadForm<3>::reduce(q, s);
    for (const Z_PrimeSymbol& symb : primes)
    {
        for (int j=3; j<=symb.power; j+=2)
        {
	  q.B_(0,2) *= symb.p;
	  q.B_(2,0) *= symb.p;
	  q.B_(1,2) *= symb.p;
	  q.B_(2,1) *= symb.p;
	  q.B_(2,2) *= symb.p * symb.p;
        }
    }
    q = Z_QuadForm<3>::reduce(q, s);

    // Do one final verification that the symbols are correct.
    Z x = q.B_(0,0) * q.B_(1,1) - q.B_(0,1) * q.B_(1,0);
    mask = 1LL << primes.size();
    for (const Z_PrimeSymbol& symb : primes)
    {
        mask >>= 1;
        int sign = (target & mask) ? -1 : 1;
        if (Z_Math::hilbert_symbol(-x, -disc, symb.p) != sign)
        {
            throw std::runtime_error("Hilbert symbols do not check out. How did this happen?");
        }
    }

    return q;
}
