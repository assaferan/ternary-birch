// implementation file for header Genus.h

template<typename R, size_t Rank>
Z Genus<R, Rank>::get_mass(const QuadForm<R, Rank>& q,
                  const std::vector<PrimeSymbol<R>>& symbols)
{
   if (Rank == 3) {
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
     // !! TODO - complete here
     return 0;
   }
}

template<typename R, size_t Rank>
Genus<R,Rank>::Genus(const QuadForm<R, Rank>& q,
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
  for (size_t n=1; n<num_conductors; n++)
    {
      if (n == 2*mask)
	{
	  ++bits;
	  mask = 1LL << bits;
	}
      R value = this->prime_divisors[bits] * this->conductors[n ^ mask];
      this->conductors.push_back(value);
    }
  
  GenusRep<R> rep;
  rep.q = q;
  rep.p = 1;
  rep.parent = -1;

  // Set the mass as a multiple of 24, as this is the largest integer
  // that can appear in its denominator. This value is used to determine
  // when the genus has been fully populated.
  this->mass_x24 = this->get_mass(q, symbols);

  // The mass provides a reasonable estimate for the size of the genus
  // since most isometry classes typically have trivial automorphism
  // group.
  Z64 estimated_size = ceil(mpz_get_d(this->mass_x24.get_mpz_t()) / 24.0);
  auto *ptr = new HashMap<GenusRep<R>>(estimated_size);
  this->hash = std::unique_ptr<HashMap<GenusRep<R>>>(ptr);
  this->hash->add(rep);

  // The spinor primes hash table, used to identify the primes used in
  // constructing the genus representatives.
  auto *ptr2 = new HashMap<W16>();
  this->spinor_primes = std::unique_ptr<HashMap<W16>>(ptr2);

  Z sum_mass_x24 = (48 / QuadForm<R>::num_automorphisms(q));

  Z p = 1;
  W16 prime = 1;
  
  // A temporary placeholder for the genus representatives before they
  // are fully built.
  GenusRep<R> foo;

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
	  const QuadForm<R>& mother = this->hash->get(current).q;
	  NeighborManager<W16,W32,R> manager(mother, GF);

#ifdef DEBUG
	  // Build the affine quadratic form for debugging purposes.
	  W16_QuadForm qp = mother.mod(GF);
#endif
	  
	  for (W16 t=0; !done && t<=prime; t++)
	    {
#ifdef DEBUG
	      // Verify that the appropriate vector is isotropic.
	      W16_Vector3 vec = manager.isotropic_vector(t);
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
	      foo.q = QuadForm<R>::reduce(foo.q, foo.s);
	      foo.p = prime;
	      foo.parent = current;
	      
	      bool added = this->hash->add(foo);
	      if (added)
		{
		  const GenusRep<R>& temp = this->hash->last();
		  sum_mass_x24 += 48 / QuadForm<R>::num_automorphisms(temp.q);
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
  for (size_t n=0; n<this->hash->size(); n++)
    {
      GenusRep<R>& rep = this->hash->at(n);
      
      // Only compute composite isometries if we are not considering the
      // mother form.
      if (n)
	{
	  GenusRep<R>& parent = this->hash->at(rep.parent);
	  
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
      const std::vector<Isometry<R>>& auts = QuadForm<R>::proper_automorphisms(rep.q);
      std::vector<bool> ignore(this->conductors.size(), false);
      for (const Isometry<R>& s : auts)
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
      
      int num = QuadForm<R>::num_automorphisms(rep.q);
      for (size_t k=0; k<num_conductors; k++)
	{
	  if (!ignore[k])
	    {
	      this->lut_positions[k][n] = this->dims[k];
	      this->num_auts[k].push_back(num);
	    }
	  this->dims[k] += (ignore[k] ? 0 : 1);
	}
    }
}
