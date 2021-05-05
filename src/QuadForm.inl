// Implementation of templated functions from QuadForm.h
#include "Math.h"

// c-tors
template<typename R, size_t n>
QuadForm_Base<R, n>::QuadForm_Base(const typename
				   QuadForm_Base<R,n>::SymVec& coeffs)
{
  size_t idx = 0;
  for (size_t row = 0; row < n; row++)
    {
      for (size_t col = 0; col < row; col++)
	{	
	  this->B_(row,col) = coeffs[idx];
	  this->B_(col,row) = coeffs[idx++];
	}
      this->B_(row,row) = coeffs[idx++];
    }
  this->is_reduced_ = false;
}

// assignment
template<typename R, size_t n>
QuadForm_Base<R,n>&
QuadForm_Base<R, n>::operator=(const QuadForm_Base<R,n> & other)
{
  if (this != &other) {
    this->B_ = other.B_;
    this->is_reduced_ = other.is_reduced_;
    this->num_aut_ = other.num_aut_;
  }
  return *this;
}

// When n is odd we return the half-discriminant
template<typename R, size_t n>
R QuadForm_Base<R, n>::discriminant(void) const
{
  R det = this->B_.determinant();
  return (n % 2 == 0) ? det : det/2;
}

// This is somewhat of a duplicate for cholesky,
// but this one keeps everything integral.
// Maybe replace cholesky in SquareMatrix with this version.
template<typename R, size_t n>
Vector<R, n> QuadForm_Base<R, n>::orthogonalize_gram() const
{
  Vector<R, n> D;
  SquareMatrix<R, n> L;
  R prod_diag = 1;
  R d, inner_sum;
  // This works but inefficiently - for some reason we get O(n^4) operations.
  // !! TODO - check it out later
  // Oh I see - we should do the L update in two passes...
  for (size_t i = 0; i < n; i++)
    {
      L(i, i) = prod_diag;
      d = prod_diag;
      for (size_t j = 0; j < i; j++)
	{
	  L(i, j) = 0;
	  for (size_t k = 0; k < i; k++)
	    {
	      inner_sum = 0;
	      for (size_t r = 0; r <= k; r++)
		inner_sum += L(k, r)*(this->B_(i,r))*L(k,j);
	      inner_sum *= -L(i, i) / D[k];
	      L(i,j) += inner_sum;
	    }
	  d = gcd(d, L(i, j));
	}
      for (size_t j = 0; j <= i; j++)
	L(i,j) /= d;
      D[i] = 0;
      for (size_t j = 0; j <= i; j++)
	for (size_t k = 0; k <= i; k++)
	  D[i] += L(i, j)*(this->B_(j,k))*L(i, k);
      prod_diag = lcm(prod_diag, D[i]);
    }

  // Recall that this is an even lattice, so all entries in D
  // are even, and we are more interested in their half values,
  // which corresponds to the quadratic form.
  for (size_t i = 0; i < n; i++)
    D[i] /= 2;
  // std::cout<< "L=" << std::endl << QuadForm(L) << std::endl;
  
  return D;
}

template<typename R, size_t n>
int QuadForm_Base<R,n>::hasse(const Vector<R, n> & D, const R & p)
{
  int hasse = 1;
  R prod = 1;
  for (size_t i = 0; i < n; i++)
    prod *= D[i];
  for (size_t i = 0; i < n-1; i++)
    {
      prod /= D[i];
      hasse *= Math<R>::hilbert_symbol(D[i], prod, p);
    }
  return hasse;
}


template<typename R, size_t n>
R QuadForm_Base<R, n>::invariants(std::set<R> & F, size_t& I) const
{
  Vector<R, n> D = this->orthogonalize_gram();
  std::set<R> P;
  F.clear();
  I = 0;
  
  P.insert(2);
  for (size_t i = 0; i < n; i++)
    {
      if (D[i] < 0) I++;
      std::vector< std::pair<R, size_t> > facs = Math<R>::factorization(D[i]);
      for (std::pair<R, size_t> fa : facs)
	  if (fa.second % 2 == 1)
	    P.insert(fa.first);
    }
  for (R p : P)
     if (Hasse(D,p) == -1) F.insert(p);

  R prod = 1;
  for (size_t i = 0; i < n; i++)
    prod *= D[i];
  
  return prod;
}

template<typename R, size_t n>
R QuadForm_Base<R, n>::invariants(std::set<std::pair<R, int> > & F, size_t& I) const
{
  Vector<R, n> D = this->orthogonalize_gram();
  std::set<R> P;
  F.clear();
  I = 0;
  
  P.insert(2);
  for (size_t i = 0; i < n; i++)
    {
      if (D[i] < 0) I++;
      std::vector< std::pair<R, size_t> > facs = Math<R>::factorization(D[i]);
      for (std::pair<R, size_t> fa : facs)
	  if (fa.second % 2 == 1)
	    P.insert(fa.first);
    }
  for (R p : P)
    F.insert(std::make_pair(p, hasse(D,p)));

  R prod = 1;
  for (size_t i = 0; i < n; i++)
    prod *= D[i];
  
  return prod;
}

template<typename R, size_t n>
typename QuadForm_Base<R, n>::jordan_data
QuadForm_Base<R, n>::jordan_decomposition(const R & p) const
{
  bool even = (p == 2);
  SquareMatrix< Rational<R>, n> S = SquareMatrix< Rational<R>, n>::identity();
  SquareMatrix< Rational<R>, n> G;
  Matrix< Rational<R> > F(n,n);
  
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      F(i,j) = this->B_(i,j);
  
  size_t k = 0;
  // virtually infinity
  size_t old_val = 0xffffffff;
  std::vector<size_t> blocks;
  jordan_data jordan;
  
  while (k < n)
    {
      // std::cerr << "k = " << k << std::endl;
      // G = SFS^t
      // !! TODO - can we write simply G = S*(this->B_)*S.transpose() ?
     for (size_t i = 0; i < n; i++)
       for (size_t j = 0; j < n; j++)
	 G(i, j) = SquareMatrix<R,n>::inner_product(this->B_, S, i, j);
     /*
     std::cerr << "G = " << std::endl;
     pretty_print<R,n>(std::cerr,G);
     */
     size_t ii = k;
     // infty
     size_t m = 0xffffffff;
     
     for (size_t i = k; i < n; i++)
       {
	 if (G(i,i) != 0) {
	   size_t val = Math<R>::valuation(G(i,i), p);
	   if (val < m)
	     {
	       m = val;
	       ii = i;
	     }
	 }
       }
     std::pair<size_t, size_t> i_pair = std::make_pair(ii, ii);
     for (size_t i = k; i < n; i++)
       for (size_t j = i+1; j < n; j++)
	 {
	   if (G(i,j) != 0) {
	     size_t tmp = Math<R>::valuation(G(i,j), p);
	     if (tmp < m)
	       {
		 m = tmp;
		 i_pair.first = i;
		 i_pair.second = j;
	       }
	   }
	 }
     /*
     std::cerr << "i_pair = (" << i_pair.first << "," << i_pair.second << ")";
     std::cerr << std::endl << "m = " << m << std::endl;
     */
     if (m != old_val)
       {
	 blocks.push_back(k);
	 old_val = m;
	 jordan.exponents.push_back(m);
       }
     /*
     std::cerr << "blocks = ";
     pretty_print<size_t>(std::cerr, blocks);
     std::cerr << "jordan.exponents = ";
     pretty_print<size_t>(std::cerr, jordan.exponents);
     */
     if ((even) && (i_pair.first != i_pair.second))
       {
	 S.swap_rows(i_pair.first, k);
	 S.swap_rows(i_pair.second, k+1);
	 
	 // T12 = S[k]*F*S[k+1]^t
	 Rational<R> T12 =
	   SquareMatrix<R,n>::inner_product(this->B_, S, k, k+1);

	 // multiply S[k] by p^val(T12,p)/T12
	 // Check whether we have to change to rational here
	 for (size_t i = 0; i < n; i++)
	   S(k,i) *= (1 << Math<R>::valuation(T12, p)) / T12;
	 Rational<R> T11 = SquareMatrix<R,n>::inner_product(this->B_, S, k, k);
	 Rational<R> T22 =
	   SquareMatrix<R,n>::inner_product(this->B_, S, k+1, k+1);
	 T12 = SquareMatrix<R,n>::inner_product(this->B_, S, k, k+1);
	 Rational<R> d = T11*T22-T12*T12;
	 for (size_t l = k+2; l < n; l++)
	   {
	     Rational<R> tl =
	       T12*SquareMatrix<R,n>::inner_product(this->B_,S,k+1,l) -
	       T22*SquareMatrix<R,n>::inner_product(this->B_,S,k,l);
	     Rational<R> ul =
	       T12*SquareMatrix<R,n>::inner_product(this->B_,S,k,l) -
	       T11*SquareMatrix<R,n>::inner_product(this->B_,S,k+1,l);
	     for (size_t i = 0; i < n; i++)
	       S(l,i) += (tl/d)*S(k,i) + (ul/d)*S(k+1,i);
	   }
	 k += 2;
       }
     else
       {
	 if (i_pair.first == i_pair.second) {
	   // std::cerr << "swapping rows" << std::endl;
	   S.swap_rows(i_pair.first, k);
	   
	   /*
	   std::cerr << "S = " << std::endl;
	   pretty_print<R,n>(std::cerr, S);
	   */
	 }
	 else
	   {
	     // std::cerr << "adding rows" << std::endl;
	     S.add_row(i_pair.first, i_pair.second, 1);
	     
	     /*
	     std::cerr << "S = " << std::endl;
	     pretty_print<R,n>(std::cerr, S);
	     std::cerr << "swapping rows" << std::endl;
	     */
	     S.swap_rows(i_pair.first, k);
	     /*
	     std::cerr << "S = " << std::endl;
	     pretty_print<R,n>(std::cerr, S);
	     */
	   }
	 Rational<R> nrm = SquareMatrix<R,n>::inner_product(this->B_, S, k, k);
	
	 // std::cerr << "nrm = " << nrm << std::endl;
	 
	 Rational<R> X[n];
	 for (size_t i = 0; i < n; i++)
	   X[i] = SquareMatrix<R,n>::inner_product(this->B_, S, k, i);
	 /*
	 std::cerr << "X = ";
	 pretty_print<Rational<R> ,n>(std::cerr, X);
	 */
	 for (size_t l = k+1; l < n; l++)
	     for (size_t i = 0; i < n; i++)
	       S(l,i) -= X[l]/nrm * S(k, i);
	 /*
         std::cerr << "S = " << std::endl;
	 pretty_print<R,n>(std::cerr, S);
	 */
	 k += 1;
       }
    }
  blocks.push_back(n);
  /*
  std::cerr << "blocks = ";
  pretty_print<size_t>(std::cerr, blocks);
  */
  for (size_t i = 0; i < blocks.size()-1; i++) {
    size_t nrows = blocks[i+1]-blocks[i];
    std::vector< Rational<R> > data(nrows*n);
    size_t idx = 0;
    for (size_t row = 0; row < nrows; row++)
      for (size_t col = 0; col < n; col++)
	data[idx++] = S(blocks[i]+row, col);
    Matrix< Rational<R> > mat(data, nrows, n);
    jordan.matrices.push_back(mat);
  }
  
  for (Matrix< Rational <R> > m  : jordan.matrices) {
    /*
    std::cerr << "m = " << m << std::endl;
    std::cerr << "F = " << F << std::endl;
    std::cerr << "m^t = " << m.transpose() << std::endl;
    */
    Rational<R> tmp_rat = m(0,0)*F(0,0);

    //    std::cerr << "tmp_rat = " << tmp_rat << std::endl;
    
    Matrix< Rational<R> > tmp = m*F;
    
    //    std::cerr << "m*F = " << tmp << std::endl;

    Matrix< Rational<R> > tmp2 = m.transpose();
    Matrix< Rational<R> > tmp3 = tmp*tmp2;

    //    std::cerr << "m*F*m^t = " << tmp3 << std::endl;
    
    jordan.grams.push_back(m*F*m.transpose());
  }
  /*
  std::cerr << "jordan.matrices = " << std::endl;
  for (size_t i = 0; i < jordan.matrices.size(); i++)
    std::cerr << std::endl << jordan.matrices[i] << std::endl;

  std::cerr << "jordan.grams = " << std::endl;
  for (size_t i = 0; i < jordan.grams.size(); i++)
    std::cerr << std::endl << jordan.grams[i] << std::endl;

  std::cerr << "jordan.exponents = ";
  for (size_t i = 0; i < jordan.exponents.size(); i++)
    std::cerr << jordan.exponents[i] << " ";
  std::cerr << std::endl;
  */
  return jordan;
}

template<typename R, size_t n>
Vector<R, n> QuadForm_Base<R,n>::voronoi_bounds(void)
{
  // !! TODO - check what the real bounds are !!
  Vector<R, n> bounds;
  for (size_t i = 0; i < n; i++)
    bounds[i] = 1;
  return bounds;
}

template<typename R, size_t n>
void
QuadForm_Base<R,n>::closest_lattice_vector(SquareMatrix<R,n> &q,
					   Isometry<R,n> & iso)
{
  // !! TODO - replace Rational by finite precision (one bit precision, maybe)
  SquareMatrix<Rational<R>, n-1> H;
  Vector<Rational<R>, n-1> v;
  Isometry<R, n> g, min_g;
  SquareMatrix<R, n> x_gram;

  for (size_t i = 0; i < n-1; i++) {
    Rational<R> scalar(1, q(i,i)); 
    for (size_t j = 0; j < n-1; j++) {
      H(i,j) = scalar*q(i,j);
    }
    v[i] = scalar*q(i, n-1);
  }
  
  Vector<Rational<R>, n-1> y = H.solve(v);
  Vector<R, n-1> voronoi = QuadForm_Base<R,n-1>::voronoi_bounds();
  Vector<R, n-1> x, x_min, x_max, x_num;
  Vector<R, n-1> x_closest;
  for (size_t i = 0; i < n-1; i++)
    x_min[i] = (y[i] - voronoi[i]).ceiling();
  for (size_t i = 0; i < n-1; i++)
    x_max[i] = (y[i] + voronoi[i]).floor();
  for (size_t i = 0; i < n-1; i++)
    x_num[i] = x_max[i] - x_min[i] + 1;
  size_t num_xs = 1;
  for (size_t i = 0; i < n-1; i++)
    num_xs *= x_num[i];
  // This should be infinity
  R min_dist = 0xffffffff;
  for (size_t x_idx = 0; x_idx < num_xs; x_idx++) {
    size_t tmp = num_xs;
    for (size_t i = 0; i < n-1; i++) {
      x[i] = x_min[i] + (tmp % x_num[i]);
      tmp /= x_num[i];
    }
    for (size_t i = 0; i < n-1; i++)
      g(i,n-1) = -x[i];
    x_gram = g.transform(q, 1);
    if (x_gram(n-1,n-1) < min_dist) {
      min_dist = x_gram(n-1,n-1);
      min_g = g;
      x_closest = x;
    }
  }
  iso *= min_g;
  q = min_g.transform(q, 1);
  return;
}
  
template<typename R, size_t n>
void QuadForm_Base<R,n>::greedy(SquareMatrix<R,n>& gram, Isometry<R,n>& s)
{
  if (n == 1) return;

  // temp isometry
  Isometry<R, n> temp;

  std::pair< R, size_t > perm_pair[n];
  Vector<size_t, n> perm;
  do {
    for (size_t i = 0; i < n; i++)
      perm_pair[i] = std::make_pair(gram(i,i), i);
    std::sort(perm_pair, perm_pair+n);
    
    for (size_t i = 0; i < n; i++)
      perm[i] = perm_pair[i].second;

    // update isometry
    s.update_perm(perm);
    temp.update_perm(perm);
    
    // update gram
    gram = temp.transform(gram, 1);
    
    // prepare arguments for recursive call
    // !! TODO - by modifying arguments to be arbitrary R**
    // we can skip this stage.
    SquareMatrix<R, n-1> subgram;
    
    for (size_t i = 0; i < n-1; i++)
      for (size_t j = 0; j < n-1; j++)
	subgram(i,j) = gram(i,j);
    
    Isometry<R,n-1> iso0;
	
    QuadForm_Base<R,n-1>::greedy(subgram, iso0);

    Isometry<R, n> iso;
    for (size_t i = 0; i < n-1; i++)
      for (size_t j = 0; j < n-1; j++)
	iso(i,j) = iso0(i,j);
    iso(n-1,n-1) = 1;

    s = s*iso;
    // !! TODO - one can use subgram to save computations
    gram = iso.transform(gram, 1);
    closest_lattice_vector(gram, iso);
    
  } while (gram(n-1,n-1) != gram(n-2,n-2));
  return;
}

template<typename R, size_t n>
std::vector< std::vector<size_t> > QuadForm_Base<R,n>::all_perms(size_t m)
{
  std::vector< std::vector<size_t> > perms;
  if (m == 1) {
    std::vector<size_t> id(1);
    id[0] = 1;
    perms.push_back(id);
    return perms;
  }
  std::vector< std::vector<size_t> > rec_perms = all_perms(m-1);
  for (std::vector<size_t> perm : rec_perms) {
    perm.push_back(m-1);
    perms.push_back(perm);
    for (size_t i = 0; i < m-1; i++) {
      std::vector<size_t> perm_new(m);
      for (size_t j = 0; j < m-1; j++)
	perm_new[j] = perm[j];
      perm_new[m-1] = perm[i];
      perm_new[i] = m-1;
      perms.push_back(perm_new);
    }
    
  }
  return perms;
}

template<typename R, size_t n>
bool
QuadForm_Base<R,n>::permutation_reduction(SquareMatrix<R, n> & qf,
					  Isometry<R,n> & isom,
					  std::set< Isometry<R, n> > & auts)
{
  bool is_reduced = true;
  std::map<R, std::vector<size_t> > stable_sets;
  Isometry<R, n> s_final;
  SquareMatrix<R, n> q0, q1;
  q0 = qf;
  
  for (size_t i = 0; i < n; i++) {
    R val = qf(i,i);
    auto search = stable_sets.find(val);
    if (search == stable_sets.end()) {
      std::vector<size_t> empty_vec;
      stable_sets[val] = empty_vec;
    }
    stable_sets[val].push_back(i);
  }
  // !! TODO - Here I go one by one, but in magma
  // Kohel tries all possible permutations (all products)
  // Could it really matter?
  typename std::map<R, std::vector<size_t> >::const_iterator iter;
  for (iter = stable_sets.begin(); iter != stable_sets.end(); iter++) {
    R key = iter->first;
    std::vector<size_t> value = iter->second;
    std::vector< std::vector<size_t> > val_perms = all_perms(value.size());
    for (std::vector<size_t> perm : val_perms) {
      Vector<size_t, n> large_perm;
      for (size_t k = 0; k < n; k++)
	large_perm[k] = k;
      for (size_t idx = 0; idx < value.size(); idx++) {
	large_perm[value[idx]] = value[perm[idx]];
      }
      Isometry<R,n> s;
      s.update_perm(large_perm);
      q1 = s.transform(qf, 1);
      if (q1 < q0) {
	q0 = q1;
	s_final = s;
	is_reduced = false;
      }
      else if (q1 == q0) {
       auts.insert(isom.inverse()*s_final.inverse()*s*isom);
      }
    }
  }
  qf = q0;
  isom = isom*s_final;
  return is_reduced;
}

template<typename R, size_t n>
bool QuadForm_Base<R,n>::sign_normalization(SquareMatrix<R, n> & qf,
					    Isometry<R,n> & isom,
					    std::set< Isometry<R, n> > & auts)
{
  bool is_reduced = true;
  Fp<R, W16> GF2(2, 0);
  std::set< VectorFp<R, W16, n > > boundary_basis;
  std::set< std::pair<size_t, size_t> > priority_set;
  
  size_t count = 0;
  for (size_t j = 0; j < n; j++)
    for (size_t k = j+1; k < n; k++) {
      Matrix< FpElement<R, W16> > w_F2(boundary_basis.size()+1, n);
      std::set< VectorFp<R, W16, n > >::const_iterator bb_ptr;
      bb_ptr = boundary_basis.begin();
      for (size_t row = 0; row < boundary_basis.size(); row++) {
	for (size_t col = 0; col < n; col++)
	  w_F2(row, col) = (*bb_ptr)[col];
	bb_ptr++;
      }
      for (size_t col = 0; col < n; col++)
	w_F2(boundary_basis.size(), col) = GF2.mod(0);
      w_F2(boundary_basis.size(), j) = GF2.mod(1);
      w_F2(boundary_basis.size(), k) = GF2.mod(1);
      if ((w_F2.rank() > count) && (qf(j,k) != 0)) {
	priority_set.insert(std::make_pair(j,k));
	VectorFp<R, W16, n > last_row(GF2);
	for (size_t col = 0; col < n; col++)
	  last_row[col] = GF2.mod(0);
	last_row[j] = GF2.mod(1);
	last_row[k] = GF2.mod(1);
	boundary_basis.insert(last_row);
	count++;
      }
    }
  std::set< VectorFp<R, W16, n > > skew_basis(GF2);
  for (std::pair<size_t, size_t> x : priority_set) {
    VectorFp<R, W16, n> vec(GF2);
    for (size_t col = 0; col < n; col++)
      vec[col] = GF2.mod(0);
    vec[x.first] = GF2.mod(1);
    vec[x.second] = GF2.mod(1);
    if (qf(x.first, x.second) < 0) {
      vec[0] = GF2.mod(1);
    }
    skew_basis.insert(vec);
  }

  Matrix< FpElement<R, W16> > w_F2(skew_basis.size(), n);
  const auto basis_ptr = skew_basis.begin();
  for (size_t row = 0; row < skew_basis.size(); row++) {
    for (size_t col = 0; col < n; col++)
      w_F2(row, col) = (*basis_ptr)[col];
    basis_ptr++;
  }
  // !! Todo - check that this is the correct kernel
  Matrix< FpElement<R, W16> > ker = w_F2.kernel();
  Isometry<R, n> s;
  is_reduced = (ker.nrows() == 0);
  for (size_t row = 0; row < ker.nrows(); row++) {
    for (size_t i = 0; i < n; i++)
      if (ker(row, i) == 1) s(i,i) = -s(i,i);
    if (s.transform(qf, 1) == qf)
      auts.insert(isom.inverse()*s*isom);
  }
  qf = s.tranform(qf, 1);
  isom *= s;
  return is_reduced;
}

// Returns the matrix obtained by a basis permutation 
// such that QF[i,i] le QF[j,j] for all i le j. }
// !! TODO - maybe use this one also to update the automorphism group
// (by permutation). Though it seems this is covered by permutation_reduction
template<typename R, size_t n>
bool QuadForm_Base<R,n>::norm_echelon(SquareMatrix<R, n> & qf,
				      Isometry<R,n> & isom)
{
  bool is_reduced = true;
  Isometry<R,n> s, u0;
  for (size_t i = 0; i < n-1; i++) {
    if (qf(i+1,i+1) < qf(i,i)) {
      s(i+1, i+1) = 0;
      s(i, i) = 0;
      s(i,i+1) = 1;
      s(i+1, i) = 1;
      qf = s.transform(qf, 1);
      u0 = s*u0;
      is_reduced = false;
    }
  }
  if (u0.a != SquareMatrix<R,n>::identity())
    is_reduced = (is_reduced) && norm_echelon(qf, isom);
  isom *= u0;
  return is_reduced;
}

template<typename R, size_t n>
bool QuadForm_Base<R,n>::neighbor_reduction(SquareMatrix<R, n> & qf,
					    Isometry<R,n> & isom,
					    std::set< Isometry<R, n> > & auts)
{
  bool is_reduced = true;
  std::vector< std::set< Vector<R, n> > > local_neighbors(1);
  Isometry<R, n> b0;
  Vector<R, n> vec;
  vec[0] = 1;
  for (size_t i = 1; i < n; i++)
    vec[i] = 0;
  local_neighbors[0].insert(vec);
  size_t num_free = 1;
  for (size_t i = 1; i < n; i++) {
    num_free *= 3;
    std::set< Vector<R, n> > free_hood;
    for (size_t x_idx = 0; x_idx < num_free; x_idx++) {
      size_t tmp = num_free;
      Vector<R, n> x;
      for (size_t j = 0; j < i; j++) {
	x[j] = (tmp % 3) - 1;
	tmp /= 3;
      }
      x[i] = 1;
      for (size_t j = i+1; j < n; j++)
	x[j] = 0;
      R norm = Vector<R,n>::inner_product(x*qf, x);
      if (norm < qf(i,i)) {
	for (size_t j = 0; j < n; j++)
	  b0(i,j) = x[j];
	qf = b0.transform(qf, 1);
	isom *= b0;
	norm_echelon(qf, isom);
	return false;
      }
      else if (norm == qf(i,i)) {
	free_hood.insert(x);
      }
    }
    local_neighbors.push_back(free_hood);
  }

  std::map< R, std::vector<size_t> > norms;
  for (size_t i = 0; i < n; i++) {
    R val = qf(i,i);
    auto search = norms.find(val);
    if (search == norms.end()) {
      std::vector<size_t> empty_vec(0);
      norms[val] = empty_vec;
    }
    norms[val].push_back(i);
  }
  typename std::map< R, std::vector<size_t> >::const_iterator iter;
  for (iter = norms.begin(); iter != norms.end(); iter++) {
    std::vector<size_t> inds = iter->second;
    std::set< Vector<R, n> > X_old;
    std::set< Vector<R, n> > X_new;
    for (size_t i : inds) {
      std::set_union(X_old.begin(), X_old.end(),
		     local_neighbors[i].begin(), local_neighbors[i].end(),
		     std::back_inserter(X_new.begin()));
      X_old = X_new;
      X_new.clear();
    }
    for (size_t i  : inds)
      local_neighbors[i] = X_old;
  }

  size_t nbs_size = 1;
  for (size_t i = 0; i < local_neighbors.size(); i++)
    nbs_size *= local_neighbors[i].size();
  std::cerr << "Original NeighborSpace size: " << nbs_size << std::endl;
  
  std::vector< std::vector< Vector<R, n> > > neighbor_space;
  for (Vector<R, n> x : local_neighbors[0]) {
    std::vector< Vector<R, n> > singleton;
    singleton.push_back(x);
    neighbor_space.push_back(singleton);
  }
  for (size_t i = 1; i < n; i++) {
    std::vector< Vector<R, n> > ns1;
    R norm = qf(i,i);
    std::vector<size_t> inds;
    for (size_t j = 0; j < i; j++)
      if (qf(j,j) == norm) inds.push_back(j);
    for (Vector<R, n> y : local_neighbors[i]) {
      for (std::vector< Vector<R, n> > c : neighbor_space) {
	bool include = true;
	for (size_t j : inds)
	  if (c[j] != y) {
	    include = false;
	    break;
	  }
	if ((include) &&
	    (abs(Vector<R,n>::inner_product(c[i-1]*qf, y)) >= abs(qf(i,i-1)))) {
	  for (Vector<R, n> cc : c)
	    ns1.push_back(cc);
	  ns1.push_back(y);
	}
	else {
	  for (size_t j = 1; j < i; j++) {
	    if (abs(Vector<R,n>::inner_product(c[j-1]*qf, c[j])) >=
		abs(qf(j,j-1))) {
	      for (Vector<R, n> cc : c)
		ns1.push_back(cc);
	      ns1.push_back(y);
	      break;
	    }
	  }
	}
      }
    }
    neighbor_space = ns1;
  }

  std::cerr << "Reduced to NeighborSpace of size:";
  std::cerr << neighbor_space.size() << std::endl;

  // !! - TODO - we can from the beginning store c as a matrix
  for (std::vector< Vector<R, n> > c : neighbor_space) {
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	b0(i,j) = c[i][j];
    if (abs(b0.determinant()) == 1) {
      SquareMatrix<R, n> q0 = b0.transform(qf, 1);
      if (q0 < qf) {
	qf = q0;
        isom *= b0;
	is_reduced = false;
	sign_normalization(qf, isom, auts);
	return false;
      }
      else if (q0 == qf) {
	auts.insert(isom.inverse()*b0*isom);
      }
    }
  }
  return is_reduced;
}

// This generates the entire automorphism group.
// We don't really need to do this.
// implement a small version of Todd-Coxeter Schreier-Sims ?!
template<typename R, size_t n>
size_t QuadForm_Base<R,n>::generate_auts(std::set< Isometry<R, n> > & auts)
{
  size_t num_aut;
  do {
    num_aut = auts.size();
    for (Isometry<R,n> s : auts) {
      for (Isometry<R,n> t : auts) {
	if (auts.find(s*t) == auts.end())
	  auts.insert(s*t);
      }
    }
    // if this condition is fullfilled we are closed under taking
    // products. Since this is a finite group, we are done.
  } while (num_aut != auts.size());
  return num_aut;
}

template<typename R, size_t n>
size_t QuadForm_Base<R,n>::num_automorphisms() const
{
  if (this->is_reduced_) return this->num_aut_;
  SquareMatrix<R, n> qf = this->B_;
  Isometry<R, n> isom;
  std::set< Isometry<R, n> > auts;
  return i_reduce(qf, isom, auts);
}

// !! - TODO - think whether we want to save this as a member.
// Right now it seems to me that most of the time we don't need it,
// so there is no need to carry it around.
template<typename R, size_t n>
std::set<Isometry<R, n>> QuadForm_Base<R,n>::proper_automorphisms() const
{
  SquareMatrix<R, n> qf = this->B_;
  Isometry<R, n> isom;
  std::set< Isometry<R, n> > auts;
  i_reduce(qf, isom, auts);
  return auts;
}

template<typename R, size_t n>
QuadForm<R,n> QuadForm_Base<R,n>::reduce(const QuadForm<R,n> & q,
					 Isometry<R,n> & isom)
{
  std::set< Isometry<R, n> > auts;
  SquareMatrix<R, n> qf = q.bilinear_form();
  size_t num_aut = i_reduce(qf, isom, auts);
  QuadForm<R,n> q_red(qf);
  q_red.num_aut_ = num_aut;
  q_red.is_reduced_ = true;
  return q_red;
}

// !! - TODO - currently this computes automorphisms always.
// Should do it only if we don't know them
template<typename R, size_t n>
size_t QuadForm_Base<R,n>::i_reduce(SquareMatrix<R, n> & qf,
				    Isometry<R,n> & isom,
				    std::set< Isometry<R, n> > & auts)
{
  greedy(qf, isom);
  bool is_reduced;
  do {
    is_reduced = true;
    is_reduced = (is_reduced) && (permutation_reduction(qf, isom, auts));
    is_reduced = (is_reduced) && (sign_normalization(qf, isom, auts));
    is_reduced = (is_reduced) && (neighbor_reduction(qf, isom, auts));
  } while (!is_reduced);
  return generate_auts(auts);
}

template<typename R, size_t n>
std::ostream& operator<<(std::ostream& os, const QuadForm<R,n>& q)
{
  os << q.bilinear_form();
  return os;
}

template<typename R, typename S, size_t n>
VectorFp<R, S ,n> QuadFormFp<R, S, n>::isotropic_vector(void) const
{
  VectorFp<R, S ,n> vec(this->GF);

  // Check the diagonal
  for (size_t i = 0; i < n; i++)
    if (this->B_(i,i) == 0) {
      vec[i] = 1;
      return vec;
    }

  if (GF->prime() == 2) return this->isotropic_vector_p2();
  
  // for now we only implement the n >= 3 case
  assert(n >= 3);

  // isometry on the submatrix of 3 first variables
  SquareMatrix<FpElement<R,S>, 3> basis =
    SquareMatrix<FpElement<R,S>, 3>::identity();
  SquareMatrix<FpElement<R, S>, 3> subM;
  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++)
      subM(i,j) = this->bilinear_form()(i,j);

  FpElement<R, S> scalar;
  
  // clear the off-diagonl entries
  for (size_t i = 0; i < 2; i++)
    for (size_t j = i+1; j < 3; j++) {
      scalar = -subM(i,j) / subM(i,i);
      subM.add_col(j, i, scalar);
      subM.add_row(j, i, scalar);
      basis.add_row(j, i, scalar);
      if (subM(j,j) == 0) {
	for (size_t k = 0; k < 3; k++)
	  vec[k] = basis(j,k);
	return vec;
      }
    }

  // Check if the first two variables alone are isotropic.
  FpElement<R,S> d = -subM(1,1)*subM(2,2);
  if (d.is_square()) {
    d = d.sqrt();
    for (size_t k = 0; k < 3; k++)
      vec[k] = basis(0,k) + (this->bilinear_form()(0,0)/d) * basis(1,k);
    return vec;
  }

  FpElement<R,S> a = subM(0,0);
  FpElement<R,S> b = subM(1,1);
  FpElement<R,S> c = subM(2,2);

  VectorFp<R, S, 2> v(GF);
  bool nonzero;
  do {
    do {
      do {
	for (size_t i = 0; i < 2; i++)
	  v[i] = GF->random();
      } while ((v[0] != 0) || (v[1] != 0));
      d = -(a*v[0]*v[0] + b*v[1]*v[1])/c;
    } while (!d.is_square());
    
    d = d.sqrt();
    nonzero = false;
    for (size_t j = 0; j < 3; j++) {
      vec[j] = v[0]*basis(0,j) + v[1]*basis(1,j) + d*basis(2,j);
      nonzero = nonzero || (vec[j] != 0);
    }
  } while (nonzero);
  return vec;
}

template<typename R, typename S, size_t n>
VectorFp<R, S, n> QuadFormFp<R, S, n>::isotropic_vector_p2(void) const
{
  VectorFp<R, S, n> vec(this->GF);
  FpElement<R, S> g;
  // If we can find a pair of orthogonal basis vectors,
  //  we can easily construct an isotropic vector.
  for (size_t i = 0; i < n-1; i++)
    for (size_t j = i+1; j < n; j++) {
      if (this->bilinear_form()(i,j) == 0) {
	g = this->bilinear_form()(j,j) / this->bilinear_form()(i,i);
	assert(g.is_square());
	g = g.sqrt();
	vec[i] = g;
	vec[j] = 1;
	return vec;
      }
    }

  // Otherwise, while the formulation is a bit more
  //  complicated, we can produce an isotropic vector
  //  by taking a linear combination of the first three
  //  basis vectors as follows:

  // Convenient references.
  FpElement<R, S> a = this->bilinear_form()(0,0);
  FpElement<R, S> b = this->bilinear_form()(1,1);
  FpElement<R, S> c = this->bilinear_form()(2,2);
  FpElement<R, S> d = this->bilinear_form()(1,2);
  FpElement<R, S> e = this->bilinear_form()(0,2);
  FpElement<R, S> f = this->bilinear_form()(0,1);

  g = (b*e*e/f/f + c + e*d/f)/a;
  assert(g.is_square());
  vec[0] = g;
  vec[1] = e/f;
  vec[2] = 1;

  return vec;
}

// general hash_value

template<size_t n>
W64 Z_QuadForm<n>::hash_value(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j <= i; j++)
      fnv = (fnv ^ mpz_get_si((this->B_(i,j)).get_mpz_t())) * FNV_PRIME;

  return fnv;
}

template<size_t n>
W64 Z64_QuadForm<n>::hash_value(void) const
{
  W64 fnv = FNV_OFFSET;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j <= i; j++)
      fnv = (fnv ^ this->B_(i,j)) * FNV_PRIME;
  return fnv;
}

