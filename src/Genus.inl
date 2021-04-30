// implementation file for header Genus.h

template<typename R, size_t dim>
std::set<R>
Genus<R, dim>::Witt_to_Hasse(const R& det,
			     const std::set<std::pair<R, int> > & finite)
{
  std::set<R> Hasse;
  int c_table[8] = {2, 1, 1, -2, -2, -1, -1, 2};
  int c_mask = c_table[dim % 8];
  R c = (c_mask / 2)*det + c_mask % 2;

  for (std::pair<R, int> x : finite)
    if (x.second != Math<R>::hilbert_symbol(-1, c, x.first))
      Hasse.insert(x.first);
  
  return Hasse;
}

template<typename R, size_t n>
Rational<Z> Genus<R, n>::local_factor(const Matrix< Rational<R> > & g,
				      const R & p)
{
  size_t m = g.ncols();
  Rational<Z> f = 1;
  Rational<Z> p_i(1,p*p);
  for (size_t i = 2; i+2 <= m; i+= 2)
    {
      f *= (1-p_i);
      p_i /= p*p;
    }
  if (m % 2 == 0) {
    f *= (1-p_i);
    return f;
  }
  size_t r = m / 2;
  R sign = (r % 2 == 0) ? 1 : -1;
  Rational<R> d = g.determinant() * sign;
  if (((Math<R>::valuation(d, p)) % 2) == 0) {
    p_i = 1;
    for (size_t i = 0; i < r; i++) p_i /= p;
    if (Math<R>::is_local_square(d, p))
      f *= 1 - p_i;
    else
      f *= 1 + p_i;
  }
  return f;
}

template<typename R, size_t n>
Rational<Z> Genus<R, n>::combine(const QuadForm<R, n>& q,
				 const R & p)
{
  assert(p != 2);
  typename QuadForm<R,n>::jordan_data jordan = q.jordan_decomposition(p);
  Rational<Z> f = 1;
  size_t e = 0;
  std::vector<size_t> ms;
  size_t m = 0;
  for (Matrix< Rational<R> > g : jordan.grams) {
    ms.push_back(g.ncols());
    m += g.ncols();
  }
  for (size_t i = 0; i < ms.size(); i++) {
    size_t t = (i == 0) ? 0 : jordan.exponents[i-1];
    e += (jordan.exponents[i]-t)*(m+1)*m/2 - jordan.exponents[i]*ms[i]*(m+1)/2;
    f *= local_factor(jordan.grams[i], p);
    m -= ms[i];
  }
  // !! We might run into trouble at 2 here
  // check if we need disc or half-disc
  // size_t v = Math<R>::valuation(q.discriminant(), p);
  Matrix<R> q_mat(q.bilinear_form());
  size_t v = Math<R>::valuation(q_mat.determinant(), p);
  if ((m % 2 == 0) && (v % 2 == 1)) {
    e += (m-1)/2;
  }
  R p_e = 1;
  for (size_t i = 0; i < e; i++) p_e *= p;
  Rational<R> denom = (1<< (jordan.grams.size()-1)) * f * p_e;
  Matrix< Rational<R> > diag =
    Matrix< Rational<R> >::diagonal_join(jordan.grams);
  return local_factor(diag, p) / denom;
}

template<typename R, size_t m>
Rational<Z> Genus<R, m>::get_mass(const QuadForm<R, m>& q,
                  const std::vector<PrimeSymbol<R>>& symbols)
{
   if (m == 3) {
        Z mass = 2 * this->disc;
        Z a = q.h() * q.h() - 4 * q.a() * q.b();
        Z b = -q.a() * this->disc;

        for (const PrimeSymbol<R>& symb : symbols)
        {
            mass *= (symb.p + Math<Z>::hilbert_symbol(a, b, symb.p));
            mass /= 2;
            mass /= symb.p;
        }

        return mass;
   }
   else {

     size_t r = m / 2;

     std::set<std::pair<R, int> > Hasse;
     // do we need the dummy? We could probably let it go
     size_t dummy;
     R det = q.invariants(Hasse, dummy);
     std::set<R> Witt = Witt_to_Hasse(det, Hasse);
     std::vector< std::pair<R, size_t> > fac = Math<R>::factorization(det);
     std::set<R> B;
     // TODO - replace these by set_union, set_difference ?
     for (std::pair<R, size_t> fa : fac)
       B.insert(fa.first);
     for (R p : Witt)
       B.insert(p);
     B.erase(2);
     for (R p : B)
       Witt.erase(p);
     size_t val2 = Math<R>::valuation(det, 2);
     
     // mass from infinity and 2
     Rational<Z> mass(1, 1<<r);    
     
     for (size_t i = 1; i < m / 2 + m % 2; i++)
       mass *= -Math<Z>::bernoulli_number(2*i)/(2*i);
     
     if (m % 2 == 1)
       {	 
         if (val2 % 2 == 1)
	   {
	     mass *= (1 << r) + ((Witt.find(2) != Witt.end()) ? -1 : 1);
	     mass /= 2;
	     Witt.erase(2);
	   }
	 if (Witt.find(2) != Witt.end())
	   { 
	     mass *= (1 << (m-1)) - 1;
	     mass /= 6;
	   }
       }
     else
       {
	 R disc = (r % 2 == 1) ? -det : det;
	 if (Math<R>::is_square(disc))
	   mass *= -Math<Z>::bernoulli_number(r)/r;
	 else
	   {
	     mass *= -Math<Z>::bernoulli_number(r, disc) / r;
	     if (r % 2 == 0)
	       mass *= -Math<Z>::bernoulli_number(r) / r;
	     if (val2 % 2 == 1)
	       {
		 mass /= 2;
		 Witt.erase(2);
	       }
	   }
	 if (Witt.find(2) != Witt.end())
	   {
	     // checking if disc is a local square at 2
	     int w = 1;
	     if ((val2 % 2 != 1) && ((disc >> val2) % 8 == 1))
		 w = -1;
	     mass *= (1<<(r-1))+w;
	     mass *= (1<<r)+w;
	     mass /= 6;
	   }
       }
     // odd places which are not unimodular or have Witt invariant -1.
     for (R p : B)
       mass *= combine(q,p);
     
     return abs(mass);
   }
}

template<typename R, size_t n>
Genus<R, n>::Genus(const QuadForm<R, n>& q,
	             const std::vector<PrimeSymbol<R>>& symbols, W64 seed)
{
  if (seed == 0)
    {
      std::random_device rd;
      seed = rd();
    }
  
  this->disc = q.discriminant();
  this->seed_ = seed;
  
  this->prime_divisors.reserve(symbols.size());
  for (const PrimeSymbol<R>& symb : symbols)
    {
      this->prime_divisors.push_back(symb.p);
    }
  
  Spinor<R> *spin = new Spinor<R>(this->prime_divisors);
  this->spinor = std::unique_ptr<Spinor<R>>(spin);
  
  if (symbols.size() > 63)
    {
      throw std::domain_error("Must have 63 or fewer prime divisors.");
    }
  
  size_t num_conductors = 1LL << symbols.size();
  
  this->conductors.reserve(num_conductors);
  this->conductors.push_back(1);
  
  size_t bits = 0;
  size_t mask = 1;
  for (size_t c=1; c<num_conductors; c++)
    {
      if (c == 2*mask)
	{
	  ++bits;
	  mask = 1LL << bits;
	}
      R value = this->prime_divisors[bits] * this->conductors[c ^ mask];
      this->conductors.push_back(value);
    }
  
  GenusRep<R, n> rep;
  rep.q = q;
  rep.p = 1;
  rep.parent = -1;

  // Set the mass as a multiple of 24, as this is the largest integer
  // that can appear in its denominator. This value is used to determine
  // when the genus has been fully populated.
  this->mass_x24 = (this->get_mass(q, symbols)).floor();

  // The mass provides a reasonable estimate for the size of the genus
  // since most isometry classes typically have trivial automorphism
  // group.
  Z64 estimated_size = ceil(mpz_get_d(this->mass_x24.get_mpz_t()) / 24.0);
  auto *ptr = new HashMap<GenusRep<R,n>>(estimated_size);
  this->hash = std::unique_ptr<HashMap<GenusRep<R,n>>>(ptr);
  this->hash->add(rep);

  // The spinor primes hash table, used to identify the primes used in
  // constructing the genus representatives.
  auto *ptr2 = new HashMap<W16>();
  this->spinor_primes = std::unique_ptr<HashMap<W16>>(ptr2);

  Z sum_mass_x24 = (48 / QuadForm<R, n>::num_automorphisms(q));

  Z p = 1;
  W16 prime = 1;
  
  // A temporary placeholder for the genus representatives before they
  // are fully built.
  GenusRep<R,n> foo;

  bool done = (sum_mass_x24 == this->mass_x24);
  while (!done)
    {
      // Get the next good prime and build the appropriate finite field.
      do
	{
	  mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
	  prime = mpz_get_ui(p.get_mpz_t());
	}
      while (this->disc % prime == 0);
      std::shared_ptr<W16_Fp> GF;
      if (prime == 2)
	GF = std::make_shared<W16_F2>(prime, this->seed_);
      else
	GF = std::make_shared<W16_Fp>(prime, this->seed_, true);
      
      size_t current = 0;
      while (!done && current < this->hash->size())
	{
	  // Get the current quadratic form and build the neighbor manager.
	  const QuadForm<R, n>& mother = this->hash->get(current).q;
	  NeighborManager<W16,W32,R,n> manager(mother, GF);

#ifdef DEBUG
	  // Build the affine quadratic form for debugging purposes.
	  W16_QuadForm<n> qp = mother.mod(GF);
#endif
	  
	  for (W16 t=0; !done && t<=prime; t++)
	    {
#ifdef DEBUG
	      // Verify that the appropriate vector is isotropic.
	      W16_Vector<n> vec = manager.isotropic_vector(t);
	      assert( qp.evaluate(vec) % prime == 0 );
#endif
	      
	      // Construct the neighbor, the isometry is stored in s.
	      foo.s.set_identity();
	      foo.q = manager.get_neighbor(t, foo.s);
	      
#ifdef DEBUG
	      // Verify neighbor discriminant matches.
	      assert( rep.q.discriminant() == mother.discriminant() );
#endif

	      // Reduce the neighbor to its Eisenstein form and add it to
	      // the hash table.
	      foo.q = QuadForm<R, n>::reduce(foo.q, foo.s);
	      foo.p = prime;
	      foo.parent = current;
	      
	      bool added = this->hash->add(foo);
	      if (added)
		{
		  const GenusRep<R,n>& temp = this->hash->last();
		  sum_mass_x24 += 48 / QuadForm<R, n>::num_automorphisms(temp.q);
		  done = (sum_mass_x24 == this->mass_x24);
		  this->spinor_primes->add(prime);
		}
	    }
	  
	  ++current;
	}
    }
  
  // Initialize the dimensions to zero, we will compute these values below.
  this->dims.resize(num_conductors, 0);
  
  // Create the lookup table values for each genus rep at each conductor.
  size_t genus_size = this->hash->size();
  this->lut_positions.resize(num_conductors, std::vector<int>(genus_size, -1));
  this->num_auts.resize(num_conductors);

  // The genus rep isometries were initialized only to contain the
  // isometry between the parent and its child, we now want to update
  // these isometries so that they are rational isometries between the
  // "mother" quadratic form and the genus rep.
  for (size_t idx=0; idx<this->hash->size(); idx++)
    {
      GenusRep<R,n>& rep = this->hash->at(idx);
      
      // Only compute composite isometries if we are not considering the
      // mother form.
      if (idx)
	{
	  GenusRep<R,n>& parent = this->hash->at(rep.parent);
	  
	  // Construct the isometries to/from the mother quadratic form.
	  rep.sinv = rep.s.inverse(rep.p);
	  rep.sinv = rep.sinv * parent.sinv;
	  rep.s = parent.s * rep.s;

	  // Copy the numerators, and increment the genus rep prime.
	  rep.es = parent.es;
	  ++rep.es[rep.p];

#ifdef DEBUG
	  R scalar = birch_util::my_pow(rep.es);
	  scalar *= scalar;
	  
	  // Verify that s is an isometry from the mother form to the rep,
	  // and that sinv is an isometry from the rep to the mother form.
	  assert( rep.s.is_isometry(q, rep.q, scalar) );
	  assert( rep.sinv.is_isometry(rep.q, q, scalar) );
#endif
	}
      
      // Determine which subspaces this representative contributes.
      const std::vector<Isometry<R, n>>& auts = QuadForm<R, n>::proper_automorphisms(rep.q);
      std::vector<bool> ignore(this->conductors.size(), false);
      for (const Isometry<R,n>& s : auts)
	{
	  Z64 vals = this->spinor->norm(rep.q, s, 1);
	  
	  for (size_t k=0; k<num_conductors; k++)
	    {
	      if (!ignore[k] && (birch_util::popcnt(vals & k) & 1))
		{
		  ignore[k] = true;
		}
	    }
	}
      
      int num = QuadForm<R, n>::num_automorphisms(rep.q);
      for (size_t k=0; k<num_conductors; k++)
	{
	  if (!ignore[k])
	    {
	      this->lut_positions[k][idx] = this->dims[k];
	      this->num_auts[k].push_back(num);
	    }
	  this->dims[k] += (ignore[k] ? 0 : 1);
	}
    }
}

template<typename R, size_t n>
template<typename T>
Genus<R, n>::Genus(const Genus<T, n>& src)
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
  this->hash = std::unique_ptr<HashMap<GenusRep<R,n>>>(new HashMap<GenusRep<R,n>>(src.hash->size()));
  for (const GenusRep<T, n>& rep : src.hash->keys())
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

template<typename R, size_t n>
Eigenvector<R> Genus<R,n>::eigenvector(const std::vector<Z32>& vec,
					  const R& conductor) const
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
  
  for (size_t idx=0; idx<fulldim; idx++)
    {
      if (lut[idx] != -1)
	{
	  temp[idx] = vec[lut[idx]];
	}
    }

  return Eigenvector<R>(std::move(temp), k);
}

template<typename R, size_t n>
std::vector<Z32>
Genus<R, n>::eigenvalues(EigenvectorManager<R, n>& vector_manager,
			    const R& p) const
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

template<typename R, size_t n>
template<typename S, typename T>
std::vector<Z32>
Genus<R,n>::_eigenvectors(EigenvectorManager<R,n>& vector_manager,
			     std::shared_ptr<Fp<S,T>> GF, const R& p) const
{
  std::vector<Z32> eigenvalues(vector_manager.size());

  S prime = GF->prime();

  const GenusRep<R,n>& mother = this->hash->get(0);

  const Z32 *stride_ptr = vector_manager.strided_eigenvectors.data();

  size_t num_indices = vector_manager.indices.size();
  for (size_t index=0; index<num_indices; index++)
    {
      size_t npos = static_cast<size_t>(vector_manager.indices[index]);
      const GenusRep<R,n>& cur = this->hash->get(npos);
      NeighborManager<S,T,R,n> neighbor_manager(cur.q, GF);
      for (W64 t=0; t<=prime; t++)
	{
	  GenusRep<R,n> foo = neighbor_manager.get_reduced_neighbor_rep((S)t);
	  
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
	      const GenusRep<R,n>& rep = this->hash->get(rpos);
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

template<typename R, size_t n>
std::map<R,std::vector<std::vector<int>>>
Genus<R, n>::hecke_matrix_sparse_internal(const R& p) const
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

  const GenusRep<R,n>& mother = this->hash->keys()[0];
  size_t num_reps = this->size();
  for (size_t idx=0; idx<num_reps; idx++)
    {
      const GenusRep<R,n>& cur = this->hash->get(idx);
      NeighborManager<W16,W32,R,n> manager(cur.q, GF);

      for (W16 t=0; t<=prime; t++)
	{
	  GenusRep<R,n> foo = manager.get_reduced_neighbor_rep(t);

#ifdef DEBUG
	  assert( foo.s.is_isometry(cur.q, foo.q, p*p) );
#endif

	  size_t r = this->hash->indexof(foo);

#ifdef DEBUG
	  assert( r < this->size() );
#endif

	  W64 spin_vals;
	  if (r == idx)
	    {
	      spin_vals = this->spinor->norm(foo.q, foo.s, p);
	    }
	  else
	    {
	      const GenusRep<R,n>& rep = this->hash->get(r);
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
	  int npos = lut[idx];
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

template<typename R, size_t n>
std::map<R,std::vector<int>>
Genus<R, n>::hecke_matrix_dense_internal(const R& p) const
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

  const GenusRep<R,n>& mother = this->hash->keys()[0];
  size_t num_reps = this->size();

  // Create hash tables for storing isotropic vectors to be skipped
  // at later iterations.
  std::vector<HashMap<W16_Vector<n> >> vector_hash(num_reps);

  for (size_t idx=0; idx<num_reps; idx++)
    {
      const GenusRep<R,n>& cur = this->hash->get(idx);
      NeighborManager<W16,W32,R,n> manager(cur.q, GF);

      for (W16 t=0; t<=prime; t++)
	{
	  GenusRep<R,n> foo;
	  W16_Vector<n> vec = manager.isotropic_vector(t);
	  vec.x = GF->mod(vec.x);
	  vec.y = GF->mod(vec.y);
	  vec.z = GF->mod(vec.z);

	  // If this vector has already been identified, skip it. This
	  // avoids unnecessarily computing the neighbor, reducing it,
	  // and testing for isometry. The Hermitian symmetry property
	  // of the Hecke matrix will account for this once we finish
	  // processing neighbors.
	  if (vector_hash[idx].exists(vec)) continue;

	  // Build the neighbor and reduce it.
	  foo.q = manager.build_neighbor(vec, foo.s);
	  foo.q = QuadForm<R, n>::reduce(foo.q, foo.s);

#ifdef DEBUG
	  assert( foo.s.is_isometry(cur.q, foo.q, p*p) );
#endif

	  size_t r = this->hash->indexof(foo);

#ifdef DEBUG
	  assert( r < this->size() );
#endif

	  W64 spin_vals;
	  if (r > idx)
	    {
	      W16_Vector<n> result = manager.transform_vector(foo, vec);
	      vector_hash[r].add(result);

	      const GenusRep<R,n>& rep = this->hash->get(r);
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
	  else if (r == idx)
	    {
	      spin_vals = this->spinor->norm(foo.q, foo.s, p);
	    }
	  else continue;

	  all_spin_vals.push_back((r << num_primes) | spin_vals);
	}

      for (size_t k=0; k<num_conductors; k++)
	{
	  const std::vector<int>& lut = this->lut_positions[k];
	  int npos = lut[idx];
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
