// Implementation of templated functions from QuadForm.h
#include "birch.h"
#include "Fp.h"
#include "Isometry.h"
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
    this->num_aut_init_ = other.num_aut_init_;
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
	  for (size_t k = j; k < i; k++)
	    {
	      inner_sum = 0;
	      for (size_t r = 0; r <= k; r++)
		inner_sum += L(k, r)*(this->B_(i,r))*L(k,j);
	      inner_sum *= -L(i, i) / D[k];
	      L(i,j) += inner_sum;
	    }
	  d = Math<R>::gcd(d, L(i, j));
	}
      for (size_t j = 0; j <= i; j++)
	L(i,j) /= d;
      D[i] = 0;
      for (size_t j = 0; j <= i; j++)
	for (size_t k = 0; k <= i; k++)
	  D[i] += L(i, j)*(this->B_(j,k))*L(i, k);
      prod_diag = Math<R>::lcm(prod_diag, D[i]);
      for (size_t j = i+1; j < n; j++)
	L(i, j) = 0;
    }

  // Recall that this is an even lattice, so all entries in D
  // are even, and we are more interested in their half values,
  // which corresponds to the quadratic form.
  for (size_t i = 0; i < n; i++)
    D[i] /= 2;
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "L=" << std::endl << L << std::endl;
#endif
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
     if (hasse(D,p) == -1) F.insert(p);

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
#ifdef DEBUG_LEVEL_FULL
      std::cerr << "k = " << k << std::endl;
#endif
      // G = SFS^t
      // !! TODO - can we write simply G = S*(this->B_)*S.transpose() ?
     for (size_t i = 0; i < n; i++)
       for (size_t j = 0; j < n; j++)
	 G(i, j) = SquareMatrix<R,n>::inner_product(this->B_, S, i, j);
#ifdef DEBUG_LEVEL_FULL
     std::cerr << "G = " << std::endl;
     G.pretty_print(std::cerr);
#endif
     size_t ii = k;
     // infty
     size_t m = 0xffffffff;

     R zero = Math<R>::zero();
     
     for (size_t i = k; i < n; i++)
       {
	 if (G(i,i) != zero) {
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
	   if (G(i,j) != zero) {
	     size_t tmp = Math<R>::valuation(G(i,j), p);
	     if (tmp < m)
	       {
		 m = tmp;
		 i_pair.first = i;
		 i_pair.second = j;
	       }
	   }
	 }
     
#ifdef DEBUG_LEVEL_FULL
     std::cerr << "i_pair = (" << i_pair.first << "," << i_pair.second << ")";
     std::cerr << std::endl << "m = " << m << std::endl;
#endif
     
     if (m != old_val)
       {
	 blocks.push_back(k);
	 old_val = m;
	 jordan.exponents.push_back(m);
       }
     
#ifdef DEBUG_LEVEL_FULL
     std::cerr << "blocks = " << blocks << std::endl;
     std::cerr << "jordan.exponents = " << jordan.exponents << std::endl;
#endif
     
     if ((even) && (i_pair.first != i_pair.second))
       {
	 S.swap_rows(i_pair.first, k);
	 S.swap_rows(i_pair.second, k+1);
	 
	 // T12 = S[k]*F*S[k+1]^t
	 Rational<R> T12 =
	   SquareMatrix<R,n>::inner_product(this->B_, S, k, k+1);

	 // multiply S[k] by p^val(T12,p)/T12
	 // Check whether we have to change to rational here
	 for (size_t i = 0; i < n; i++) {
	   R val = (1 << Math<R>::valuation(T12, p));
	   S(k,i) *= val / T12;
	 }
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
	   
#ifdef DEBUG_LEVEL_FULL
	   std::cerr << "swapping rows" << std::endl;
#endif
	   S.swap_rows(i_pair.first, k);
	   
#ifdef DEBUG_LEVEL_FULL
	   std::cerr << "S = " << std::endl;
	   S.pretty_print(std::cerr);
#endif
	 }
	 else
	   {
	     R one = Math<R>::one();
	     // std::cerr << "adding rows" << std::endl;
	     S.add_row(i_pair.first, i_pair.second, one);
	     
#ifdef DEBUG_LEVEL_FULL
	     std::cerr << "S = " << std::endl;
	     S.pretty_print(std::cerr);
	     std::cerr << "swapping rows" << std::endl;
#endif
	     S.swap_rows(i_pair.first, k);
#ifdef DEBUG_LEVEL_FULL
	     std::cerr << "S = " << std::endl;
	     S.pretty_print(std::cerr);
#endif
	   }
	 Rational<R> nrm = SquareMatrix<R,n>::inner_product(this->B_, S, k, k);

#ifdef DEBUG_LEVEL_FULL
	 std::cerr << "nrm = " << nrm << std::endl;
#endif
	 
	 Rational<R> X[n];
	 for (size_t i = 0; i < n; i++)
	   X[i] = SquareMatrix<R,n>::inner_product(this->B_, S, k, i);
	 
#ifdef DEBUG_LEVEL_FULL
	 std::cerr << "X = " << X << std::endl;;
#endif
	 for (size_t l = k+1; l < n; l++)
	     for (size_t i = 0; i < n; i++)
	       S(l,i) -= X[l]/nrm * S(k, i);
	 
#ifdef DEBUG_LEVEL_FULL
         std::cerr << "S = " << std::endl;
	 S.pretty_print(std::cerr);
#endif
	 k += 1;
       }
    }
  blocks.push_back(n);
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "blocks = " << blocks << std::endl;
#endif
  
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
#ifdef DEBUG_LEVEL_FULL
    std::cerr << "m = " << m << std::endl;
    std::cerr << "F = " << F << std::endl;
    std::cerr << "m^t = " << m.transpose() << std::endl;
    Rational<R> tmp_rat = m(0,0)*F(0,0);
    
    std::cerr << "tmp_rat = " << tmp_rat << std::endl;
    Matrix< Rational<R> > tmp = m*F;
    
    std::cerr << "m*F = " << tmp << std::endl;
    Matrix< Rational<R> > tmp2 = m.transpose();
    Matrix< Rational<R> > tmp3 = tmp*tmp2;

    std::cerr << "m*F*m^t = " << tmp3 << std::endl;
#endif
    jordan.grams.push_back(m*F*m.transpose());
  }
#ifdef DEBUG_LEVEL_FULL
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
#endif
  return jordan;
}

template<typename R, size_t n>
Vector<R, n-1> QuadForm_Base<R,n>::voronoi_bounds(size_t dim)
{
  // !! TODO - check what the real bounds are !!
  Vector<R, n-1> bounds;
  for (size_t i = 0; i < dim; i++)
    bounds[i] = 1;
  return bounds;
}

template<typename R, size_t n>
void
QuadForm_Base<R,n>::closest_lattice_vector(SquareMatrix<R,n> &q,
					   Isometry<R,n> & iso,
					   size_t dim)
{
  // !! TODO - replace Rational by finite precision (one bit precision, maybe)
  SquareMatrix<Rational<R>, n-1> H = SquareMatrix<Rational<R>, n-1>::identity();
  Vector<Rational<R>, n-1> v;
  Isometry<R, n> g, min_g;
  SquareMatrix<R, n> x_gram;
  SquareMatrix<R, n-1> H_int;
  Vector<R, n-1> v_int;

#ifdef DEBUG
  std::cerr << "finding closest_lattice_vector with gram:" << std::endl;
  q.pretty_print(std::cerr, dim);
#endif
  /*
  R max_v = q(0,0);
  size_t max_i = 0;
  for (size_t i = 1; i < dim-1; i++) {
    if (max_v * q(i,i) < q(i, dim-1) * q(max_i, max_i)) {
      max_i = i;
      max_v = q(i, dim-1);
    }
  }
  
  int r = Math<R>::log2(max_v)-Math<R>::log2(q(max_i,max_i));
  if (max_v > (1 << r)*q(max_i, max_i)) r++;
  */
  for (size_t i = 0; i < dim-1; i++) {
    Rational<R> scalar(1, q(i,i));
    for (size_t j = 0; j < dim-1; j++) {
      H(i,j) = scalar*q(i,j);
      // H_int(i,j) = (q(i,j) << r) / q(i,i);
      H_int(i,j) = q(i,j);
    }
    v[i] = scalar*q(i,dim-1);
    //  v_int[i] = (q(i, dim-1) << (r+1)) / q(i,i);
    v_int[i] = q(i,dim-1);
  }

  H_int = H_int.adjugate(dim-1);
  
#ifdef DEBUG
  std::cerr << "H = " << std::endl;
  H.pretty_print(std::cerr, dim-1);

  std::cerr << "v = " << std::endl;
  v.pretty_print(std::cerr, dim-1);

  std::cerr << "H_int = " << std::endl;
  H_int.pretty_print(std::cerr, dim-1);

  std::cerr << "v_int = " << std::endl;
  v_int.pretty_print(std::cerr, dim-1);
#endif
  
  Vector<Rational<R>, n-1> y = v*H.inverse().transpose();

#ifdef DEBUG
  std::cerr << "y = " << std::endl;
  y.pretty_print(std::cerr, dim-1);

  // need to compute the pseudo-inverse here
  // y_int is 2*y mod 1 (so we have one-bit accuracy)
  
  Vector<R, n-1> y_int = v_int*H_int.transpose();
  
  std::cerr << "y_int = " << std::endl;
  y_int.pretty_print(std::cerr, dim-1);
  
#endif
  
  Vector<R, n-1> voronoi = voronoi_bounds(dim-1);
  Vector<R, n-1> x, x_min, x_max, x_num;
  Vector<R, n-1> x_closest;
  for (size_t i = 0; i < dim-1; i++)
    x_min[i] = (y[i] - voronoi[i]).ceiling();
  for (size_t i = 0; i < dim-1; i++)
    x_max[i] = (y[i] + voronoi[i]).floor();

#ifdef DEBUG
  Vector<R, n-1> x_min_int, x_max_int;
  
  // This can be calculated more efficiently
  R det = Math<R>::zero();
  for (size_t i = 0; i < dim-1; i++)
    det += H_int(0,i)*q(i,0);
  det = (det > 0) ?  det : -det;
  for (size_t i = 0; i < dim-1; i++) {
    R tmp =  y_int[i] - det*voronoi[i];
    x_min_int[i] = ((tmp >= 0) ? tmp+det-1 : tmp)/det;
  }
  for (size_t i = 0; i < dim-1; i++) {
    R tmp =  y_int[i] + det*voronoi[i];
    x_max_int[i] = ((tmp >= 0) ? tmp : tmp-det+1)/det;
  }

  for (size_t i = 0; i < dim-1; i++) {
    assert(x_min[i] == x_min_int[i]);
    assert(x_max[i] == x_max_int[i]);
  }
#endif
  for (size_t i = 0; i < dim-1; i++)
    x_num[i] = x_max[i] - x_min[i] + 1;
  R num_xs = 1;
  for (size_t i = 0; i < dim-1; i++)
    num_xs *= x_num[i];
  // This should be infinity
  R min_dist = 0xffffffff;
  for (R x_idx = 0; x_idx < num_xs; x_idx++) {
    R tmp = x_idx;
    for (size_t i = 0; i < dim-1; i++) {
      size_t j = dim-2-i;
      x[j] = x_min[j] + (tmp % x_num[j]);
      tmp /= x_num[j];
    }
    for (size_t i = 0; i < dim-1; i++)
      g(i,dim-1) = -x[i];
    x_gram = g.transform(q);
    if (x_gram(dim-1,dim-1) < min_dist) {
      min_dist = x_gram(dim-1,dim-1);
      min_g = g;
      x_closest = x;
    }
  }
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "x_closest = " << std::endl;
  x_closest.pretty_print(std::cerr, dim-1);
#endif
  
  iso = iso*min_g;
  q = min_g.transform(q);

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "returning isometry: " << std::endl;
  iso.a.pretty_print(std::cerr, dim);
  std::cerr << "transformed gram to: " << std::endl;
  q.pretty_print(std::cerr, dim);
#endif
  return;
}

// to avoid recursive template instantiation,
// we supply a parameter defining the level of recursion
// and use only this part of the matrices
// All containers will have size n, but we will only use dim entries
template<typename R, size_t n>
void QuadForm_Base<R,n>::greedy(SquareMatrix<R,n>& gram,
				Isometry<R,n>& s,
				size_t dim)
{

#ifdef DEBUG
  Isometry<R,n> s0 = s;
  SquareMatrix<R, n> q0 = gram;
#endif
  
  if (dim == 1) return;

  // temp isometry
  Isometry<R, n> temp;

  std::pair< R, size_t > perm_pair[n];
  Vector<size_t, n> perm;
  do {
    for (size_t i = 0; i < dim; i++)
      perm_pair[i] = std::make_pair(gram(i,i), i);
    std::sort(perm_pair, perm_pair+dim);
    
    for (size_t i = 0; i < dim; i++)
      perm[i] = perm_pair[i].second;

    // this is to make sure we do not touch these rows
    for (size_t i = dim; i < n; i++)
      perm[i] = i;

    temp.set_identity();
    temp.update_perm(perm);

    // update isometry
    // s.update_perm(perm);
    s = s*temp;
    
    // update gram
    gram = temp.transform(gram);
    
#ifdef DEBUG
    assert((s0.inverse()*s).transform(q0) == gram);
#endif

    // !! - TODO - do we really need iso here
    // or could we simply pass s?
    Isometry<R,n> iso;
	
    greedy(gram, iso, dim-1);

    s = s*iso;
    // !! TODO - one can use subgram to save computations
    // This transformation already happens inside greedy(dim-1)
    //     gram = iso.transform(gram);

#ifdef DEBUG
    assert((s0.inverse()*s).transform(q0) == gram);
#endif
    
    closest_lattice_vector(gram, s, dim);

#ifdef DEBUG
    assert((s0.inverse()*s).transform(q0) == gram);
#endif
    
  } while (gram(dim-1,dim-1) < gram(dim-2,dim-2));
  return;

}

template<typename R, size_t n>
std::vector< std::vector<size_t> >
QuadForm_Base<R,n>::all_perms(size_t m)
{
  std::vector< std::vector<size_t> > perms;
  if (m == 1) {
    std::vector<size_t> id(1);
    id[0] = 0;
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
    //    R key = iter->first;
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
      q1 = s.transform(qf);
      if (q1 < q0) {
	q0 = q1;
	s_final = s;
	is_reduced = false;
      }
      else if (q1 == q0) {
	auts.insert(isom*s*s_final.inverse()*isom.inverse());
      }
    }
  }
  qf = q0;
  isom = isom*s_final;
  return is_reduced;
}

template<typename R, size_t n>
bool
QuadForm_Base<R,n>::sign_normalization_slow(SquareMatrix<R, n> & qf,
					    Isometry<R,n> & isom,
					    std::set< Isometry<R, n> > & auts)
{
  bool is_reduced = true;
  W16 prime = 2;
  std::random_device rd;
  W64 seed = rd();
  std::shared_ptr<W16_F2> GF2 = std::make_shared<W16_F2>(prime,seed);
  // !! - TODO - we don't really need set here, can use vector
  std::set< W16_VectorFp<n> > boundary_basis;
  std::set< std::pair<size_t, size_t> > priority_set;
  
  size_t count = 0;
  for (size_t j = 1; j < n; j++)
    for (size_t k = 0; k < n-j; k++) {
      W16_MatrixFp w_F2(GF2, boundary_basis.size()+1, n);
      typename std::set< W16_VectorFp<n> >::const_iterator bb_ptr;
      bb_ptr = boundary_basis.begin();
      for (size_t row = 0; row < boundary_basis.size(); row++) {
	for (size_t col = 0; col < n; col++)
	  w_F2(row, col) = (*bb_ptr)[col];
	bb_ptr++;
      }
      for (size_t col = 0; col < n; col++)
	w_F2(boundary_basis.size(), col) = GF2->mod(0);
      w_F2(boundary_basis.size(), k) = GF2->mod(1);
      w_F2(boundary_basis.size(), k+j) = GF2->mod(1);
      if ((w_F2.rank() > count) && (qf(k,k+j) != 0)) {
	priority_set.insert(std::make_pair(k,k+j));
	W16_VectorFp<n> last_row(GF2);
	for (size_t col = 0; col < n; col++)
	  last_row[col] = GF2->mod(0);
	last_row[k] = GF2->mod(1);
	last_row[k+j] = GF2->mod(1);
	boundary_basis.insert(last_row);
	count++;
      }
    }

  W16_MatrixFp w_F2(GF2, priority_set.size(), n+1);
  std::set< std::pair<size_t, size_t> >::const_iterator ps_ptr;
  ps_ptr = priority_set.begin();
  for (size_t row = 0; row < priority_set.size(); row++) {
    for (size_t col = 0; col <= n; col++)
	w_F2(row, col) = GF2->mod(0);
    
    w_F2(row, ps_ptr->first) = GF2->mod(1);
    w_F2(row, ps_ptr->second) = GF2->mod(1);

    // the affine coordinate
    if (qf(ps_ptr->first, ps_ptr->second) < 0)
      w_F2(row, n) = GF2->mod(1);
    
    ps_ptr++;
  }

  W16_MatrixFp ker = w_F2.kernel();
  // The last row of ker should now be a solution to the affine equation
  // The rows above are the kernel
#ifdef DEBUG
  for (size_t row = 0; row + 1 < ker.nrows(); row++)
    assert(ker(row, n) == GF2->mod(0));
  if (ker.nrows() >= 1)
    assert(ker(ker.nrows()-1, n) == GF2->mod(1));
#endif
  //  W16_MatrixFp tmp(GF2, ker.nrows(), ker.nrows());
  // W16_MatrixFp::row_echelon(ker, tmp);
  Isometry<R, n> s;
  is_reduced = true;
  for (size_t row = 0; row + 1 < ker.nrows(); row++) {
    is_reduced = false;
    for (size_t i = 0; i < n; i++)
      s(i,i) = (ker(row, i) + ker(ker.nrows()-1, i) == 1) ? -1 : 1;
    if (s.transform(qf) == qf) {
      auts.insert(isom*s*isom.inverse());
      is_reduced = true;
      // to be compatible with magma implementation for debugging
      for (size_t i = 0; i < n; i++) s(i,i) = 1;
    }
  }
  qf = s.transform(qf);
  isom = isom*s;
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
#ifdef DEBUG
  SquareMatrix<R,n> qf_orig = qf;
#endif
  bool is_reduced = true;
  Isometry<R,n> s, u0;
  for (size_t i = 0; i < n-1; i++) {
    if (qf(i+1,i+1) < qf(i,i)) {
      s.set_identity();
      s(i+1, i+1) = 0;
      s(i, i) = 0;
      s(i,i+1) = 1;
      s(i+1, i) = 1;
      qf = s.transform(qf);
      u0 = u0*s;
#ifdef DEBUG
      assert(u0.transform(qf_orig) == qf);
#endif
      is_reduced = false;
    }
  }
  if (u0.a != SquareMatrix<R,n>::identity())
    is_reduced = (is_reduced) && norm_echelon(qf, isom);
  isom = isom*u0;
  return is_reduced;
}

template<typename R, size_t n>
bool QuadForm_Base<R,n>::neighbor_reduction(SquareMatrix<R, n> & qf,
					    Isometry<R,n> & isom,
					    std::set< Isometry<R, n> > & auts,
					    bool calc_aut)
{
#ifdef DEBUG_LEVEL_FULL
  SquareMatrix<R,n> qf_orig = qf;
  Isometry<R,n> isom_orig = isom;
#endif
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
      size_t tmp = x_idx;
      Vector<R, n> x;
      for (size_t j = 0; j < i; j++) {
	// we separate because tmp is unsigned, which might lead to overflow
	x[j] = (tmp % 3);
	x[j]--;
	tmp /= 3;
      }
      x[i] = 1;
      for (size_t j = i+1; j < n; j++)
	x[j] = 0;
      R norm = Vector<R,n>::inner_product(x*qf, x);
      if (norm < qf(i,i)) {
	b0.set_identity();
	for (size_t j = 0; j < n; j++)
	  b0(j,i) = x[j];
	qf = b0.transform(qf);
	isom = isom*b0;
#ifdef DEBUG_LEVEL_FULL
	assert((isom_orig.inverse() * isom).transform(qf_orig) == qf);
#endif
	norm_echelon(qf, isom);
#ifdef DEBUG_LEVEL_FULL
	assert((isom_orig.inverse() * isom).transform(qf_orig) == qf);
#endif
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
  // !! TODO - here there's duplication.
  // we can simply call local_neighbors[norms[i]],
  // changin local_neighbors to depend on the norm
  for (iter = norms.begin(); iter != norms.end(); iter++) {
    std::vector<size_t> inds = iter->second;
    std::set< Vector<R, n> > X;
    for (size_t i : inds) {
      X.insert(local_neighbors[i].begin(), local_neighbors[i].end());
    }
    for (size_t i  : inds)
      local_neighbors[i] = X;
  }

  size_t nbs_size = 1;
  for (size_t i = 0; i < local_neighbors.size(); i++)
    nbs_size *= local_neighbors[i].size();
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "Original NeighborSpace size: " << nbs_size << std::endl;
#endif
  
  std::vector< std::vector< Vector<R, n> > > neighbor_space;
  for (Vector<R, n> x : local_neighbors[0]) {
    std::vector< Vector<R, n> > singleton;
    singleton.push_back(x);
    neighbor_space.push_back(singleton);
  }
  for (size_t i = 1; i < n; i++) {
    std::vector< std::vector< Vector<R, n> > > ns1;
    R norm = qf(i,i);
    std::vector<size_t> inds;
    for (size_t j = 0; j < i; j++)
      if (qf(j,j) == norm) inds.push_back(j);
    for (Vector<R, n> y : local_neighbors[i]) {
      for (std::vector< Vector<R, n> > c : neighbor_space) {
	bool include = true;
	for (size_t j : inds)
	  if (c[j] == y) {
	    include = false;
	    break;
	  }
	if ((include) &&
	    (abs(Vector<R,n>::inner_product(c[i-1]*qf, y)) >= abs(qf(i,i-1)))) {
	  c.push_back(y);
	  ns1.push_back(c);
	}
	else {
	  for (size_t j = 1; j < i; j++) {
	    if (abs(Vector<R,n>::inner_product(c[j-1]*qf, c[j])) >
		abs(qf(j,j-1))) {
	      c.push_back(y);
	      ns1.push_back(c);
	      break;
	    }
	  }
	}
      }
    }
    neighbor_space = ns1;
  }
  
#ifdef DEBUG_LEVEL_FULL
  std::cerr << "Reduced to NeighborSpace of size:";
  std::cerr << neighbor_space.size() << std::endl;
#endif
  
  // !! - TODO - we can from the beginning store c as a matrix
  for (std::vector< Vector<R, n> > c : neighbor_space) {
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	// note that we transpose
	b0(i,j) = c[j][i];
    if (abs(b0.a.determinant()) == 1) {
      SquareMatrix<R, n> q0 = b0.transform(qf);
      Isometry<R, n> u;
      std::set< Isometry<R,n> > tmp_auts;
      sign_normalization(q0, u, tmp_auts, calc_aut);
      if (q0 < qf) {
	qf = q0;
        isom = isom*b0*u;
	is_reduced = false;

#ifdef DEBUG_LEVEL_FULL
	assert((isom_orig.inverse() * isom).transform(qf_orig) == qf);
#endif
	return false;
      }
      else if ((calc_aut) && (q0 == qf)) {
	auts.insert(isom*b0*u*isom.inverse());
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
    typename std::set< Isometry<R, n> >::const_iterator iter1, iter2;
    for (iter1 = auts.begin(); iter1 != auts.end(); iter1++) {
      for (iter2 = auts.begin(); iter2 != auts.end(); iter2++) {
	Isometry<R, n> prod = (*iter1)*(*iter2);
	if (auts.find(prod) == auts.end())
	  auts.insert(prod);
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
  if (this->num_aut_init_) return this->num_aut_;
  SquareMatrix<R, n> qf = this->B_;
  Isometry<R, n> isom;
  std::set< Isometry<R, n> > auts;
  return i_reduce(qf, isom, auts, true);
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
  i_reduce(qf, isom, auts, true);
  return auts;
}

template<typename R, size_t n>
QuadForm<R,n> QuadForm_Base<R,n>::reduce(const QuadForm<R,n> & q,
					 Isometry<R,n> & isom,
					 bool calc_aut)
{
#ifdef DEBUG
  assert(q.bilinear_form().is_positive_definite());
#endif
  std::set< Isometry<R, n> > auts;
  SquareMatrix<R, n> qf = q.bilinear_form();
  size_t num_aut = i_reduce(qf, isom, auts, calc_aut);
  QuadForm<R,n> q_red(qf);
  if (calc_aut) {
    q_red.num_aut_ = num_aut;
    q_red.num_aut_init_ = true;
  }
  q_red.is_reduced_ = true;
  return q_red;
}

// !! - TODO - currently this computes automorphisms always.
// Should do it only if we don't know them
template<typename R, size_t n>
size_t QuadForm_Base<R,n>::i_reduce(SquareMatrix<R, n> & qf,
				    Isometry<R,n> & isom,
				    std::set< Isometry<R, n> > & auts,
				    bool calc_aut)
{
#ifdef DEBUG
  SquareMatrix<R, n> q0 = qf;
  Isometry<R,n> s0 = isom;
#endif
  greedy(qf, isom);
#ifdef DEBUG
  assert((s0.inverse()*isom).transform(q0) == qf);
#endif
  
  bool is_reduced;
  do {
    is_reduced = true;
    is_reduced = (permutation_reduction(qf, isom, auts)) && (is_reduced);
#ifdef DEBUG
    assert((s0.inverse()*isom).transform(q0) == qf);
    for (Isometry<R, n> s : auts) {
      assert((s0.inverse()*s*s0).transform(q0) == q0);
    }
#endif    
    is_reduced = (sign_normalization(qf, isom, auts, calc_aut)) && (is_reduced);
#ifdef DEBUG
    assert((s0.inverse()*isom).transform(q0) == qf);
    for (Isometry<R, n> s : auts) {
      assert((s0.inverse()*s*s0).transform(q0) == q0);
    }
#endif
    is_reduced = (neighbor_reduction(qf, isom, auts, calc_aut)) && (is_reduced);
#ifdef DEBUG
    assert((s0.inverse()*isom).transform(q0) == qf);
    for (Isometry<R, n> s : auts) {
      assert((s0.inverse()*s*s0).transform(q0) == q0);
    }
#endif
  } while (!is_reduced);
  if (!calc_aut)
    return auts.size();
  return generate_auts(auts);
}

template<typename R, size_t n>
std::ostream& operator<<(std::ostream& os, const QuadForm<R,n> & q)
{
  os << q.bilinear_form();
  return os;
}

template<typename R, typename S, size_t n>
bool
QuadFormFp<R, S, n>::isotropic_vector(VectorFp<R, S ,n> & vec,
				      size_t start,
				      bool deterministic) const
{
  // Check the diagonal
  for (size_t i = start; i < n; i++)
    if (this->B_Fp(i,i) == 0) {
      vec[i] = 1;
      return true;
    }

  size_t dim = n - start;

  // if the one vector was isotropic, we would have found it above.
  if (dim == 1)
    return false;

  if (GF->prime() == 2) return this->isotropic_vector_p2(vec, start);

  if (dim == 2) {
    FpElement<R, S> a = this->B_Fp(start,start);
    FpElement<R, S> b = this->B_Fp(start,start+1);
    FpElement<R, S> c = this->B_Fp(start+1,start+1);
    // The form is isotropic if and only if b^2-ac is a square.
    FpElement<R, S> d = b*b-a*c;
    // If not a square, this form is anisotropic.
    if (!d.is_square()) {
      return false;
    }
    // Since a ne 0 and the form is isotropic, we're done.
    d = d.sqrt();
    vec[start] = -((b+d)/a);
    vec[start+1] = 1;
    return true;
  }

 #ifdef DEBUG
  assert(dim >= 3);
#endif

  // isometry on the submatrix of 3 first variables
  SquareMatrixFp<R, S, 3> basis(this->GF, SquareMatrixFp<R, S, 3>::identity());
  SquareMatrixFp<R, S, 3> subM(this->GF);
  for (size_t i = 0; i < 3; i++)
    for (size_t j = 0; j < 3; j++)
      subM(i,j) = this->bilinear_form()(start+i,start+j);

  FpElement<R, S> scalar;
  
  // clear the off-diagonal entries
  for (size_t i = 0; i < 2; i++)
    for (size_t j = i+1; j < 3; j++) {
      scalar = -subM(i,j) / subM(i,i);
      subM.add_col(j, i, scalar);
      subM.add_row(j, i, scalar);
      basis.add_row(j, i, scalar);
      if (subM(j,j) == 0) {
	for (size_t k = 0; k < 3; k++)
	  vec[start+k] = basis(j,k);
	return true;
      }
    }

  // Check if the first two variables alone are isotropic.
  FpElement<R,S> d = -subM(0,0)*subM(1,1);
  if (d.is_square()) {
    d = d.sqrt();
    for (size_t k = 0; k < 3; k++)
      vec[start+k] = basis(0,k) +
	(this->bilinear_form()(start,start)/d) * basis(1,k);
    return true;
  }

  if (deterministic) {
    R p = GF->prime();
    // The quadratic form over three variables.
    QuadFormFp<R, S, 3> Q(subM);
    VectorFp<R, S, 3> v(GF);
    for (R x = 0; x < p; x++)
      for (R y = 0; y < p; y++) {
	v[0] = x;
	v[1] = y;
	v[2] = 1;
	if (Q.evaluate(v) == 0) {
	  // Found an isotropic vector, return it.
	  for (size_t j = 0; j < 3; j++) {
	    vec[start+j] = v[0]*basis(0,j) + v[1]*basis(1,j) + basis(2,j);
	  }
	  return true;
	}
      }
  }

  // If we're fine with a probabilitistic means of finding
  //  isotropic vectors, we can find them much faster.
  
  FpElement<R,S> a = subM(0,0);
  FpElement<R,S> b = subM(1,1);
  FpElement<R,S> c = subM(2,2);

  VectorFp<R, S, 2> v(this->GF);
  bool nonzero;
  do {
    do {
      do {
	for (size_t i = 0; i < 2; i++)
	  v[i] = GF->random();
      } while ((v[0] == 0) && (v[1] == 0));
      d = -(a*v[0]*v[0] + b*v[1]*v[1])/c;
    } while (!d.is_square());
    
    d = d.sqrt();
    nonzero = false;
    for (size_t j = 0; j < 3; j++) {
      vec[start+j] = v[0]*basis(0,j) + v[1]*basis(1,j) + d*basis(2,j);
      nonzero = nonzero || (vec[start+j] != 0);
    }
  } while (!nonzero);
  return true;
}

template<typename R, typename S, size_t n>
bool QuadFormFp<R, S, n>::isotropic_vector_p2(VectorFp<R, S, n> & vec,
					      size_t start) const
{
  FpElement<R, S> g;
  // If we can find a pair of orthogonal basis vectors,
  //  we can easily construct an isotropic vector.
  for (size_t i = start; i < n-1; i++)
    for (size_t j = i+1; j < n; j++) {
      if (this->bilinear_form()(i,j) == 0) {
	g = this->bilinear_form()(j,j) / this->bilinear_form()(i,i);
	assert(g.is_square());
	g = g.sqrt();
	vec[i] = g;
	vec[j] = 1;
	return true;
      }
    }

  size_t dim = n - start;
  if (dim == 1) {
    if (this->bilinear_form()(start,start) == 0) {
      vec[start] = 1;
      return true;
    }
    return false;
  }
  
  if (dim == 2) {
    for (size_t i = 0; i < 2; i++)
      if (this->bilinear_form()(start+i,start+i) == 0) {
	vec[start+i] = 1;
	return true;
      }
    FpElement<R, S> a = this->B_Fp(start,start);
    FpElement<R, S> b = this->B_Fp(start,start+1);
    FpElement<R, S> c = this->B_Fp(start+1,start+1);
    if (b == 0) {
      vec[start] = 1;
      vec[start+1] = 1;
      return true;
    }
    // In this case a = b = c = 1, so the form is anisotropic
    return false;
  }
  assert (dim >= 3);
  // Otherwise, while the formulation is a bit more
  //  complicated, we can produce an isotropic vector
  //  by taking a linear combination of the first three
  //  basis vectors as follows:

  // Convenient references.
  FpElement<R, S> a = this->bilinear_form()(start,start);
  FpElement<R, S> b = this->bilinear_form()(start+1,start+1);
  FpElement<R, S> c = this->bilinear_form()(start+2,start+2);
  FpElement<R, S> d = this->bilinear_form()(start+1,start+2);
  FpElement<R, S> e = this->bilinear_form()(start,start+2);
  FpElement<R, S> f = this->bilinear_form()(start,start+1);

  g = (b*e*e/f/f + c + e*d/f)/a;
  assert(g.is_square());
  vec[start] = g;
  vec[start+1] = e/f;
  vec[start+2] = 1;

  return true;
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

template<typename R, typename S, size_t n>
void QuadFormFp<R, S, n>::split_hyperbolic_plane(const VectorFp<R, S, n>& vec,
						 SquareMatrixFp<R, S, n>& gram,
						 SquareMatrixFp<R, S, n> & basis,
						 size_t start) const
{
  // !! TODO - see if we can get rid of this.
  // SquareMatrixFp<R, S, n> orig_basis = basis;
  // The change of basis which preserving the isometry.
  basis.set_identity();
  // Make a copy of the Gram matrix.
  gram = this->bilinear_form();
  R p = GF->prime();
  // Set the diagonal entries to zero when in characteristic 2.
  // This is because we are decomposing the associated bilinear form
  if (p == 2) {
    for (size_t i = start; i < n; i++) gram(i,i) = 0;
  }
  SquareMatrixFp<R, S, n> original_gram = gram;
  // Find the pivot of the specified vector.
  size_t pivot = start;
  while (vec[pivot] == 0) pivot++;

#ifdef DEBUG
  assert(pivot < n);
#endif 
  
  // Perform the necessary basis changes so that vec becomes the first
  //  basis vector.
  basis.multiply_row(pivot, vec[pivot]);
  gram.multiply_col(pivot, vec[pivot]);
  gram.multiply_row(pivot, vec[pivot]);
  for (size_t i = pivot+1; i < n; i++) {
    basis.add_row(pivot, i, vec[i]);
    gram.add_col(pivot, i, vec[i]);
    gram.add_row(pivot, i, vec[i]);
  }
  basis.swap_rows(start, pivot);
  gram.swap_cols(start, pivot);
  gram.swap_rows(start, pivot);

  bool is_in_radical = true;
  
  // Find a basis vector which is not orthogonal to our isotropic vector.
  size_t idx = start;
  for (; idx < n; idx++)
    if (gram(start, idx) != 0) {
      is_in_radical = false;
      break;
    }

  // If the first row is entirely zero, then this vector belongs to the
  //  radical of the form.
  
  if (is_in_radical) {
    if (p == 2) {
      // Recover the quadratic form along the diagonal.
      for (size_t i = start; i < n; i++)
	gram(i, i) = this->evaluate(basis[i]);
    }
    return;
  }
  
  // Swap this basis vector with the second basis.
  basis.swap_rows(start+1, idx);
  gram.swap_cols(start+1, idx);
  gram.swap_rows(start+1, idx);

  // Normalize the second basis vector so the (0,1)-entry is 1.
  FpElement<R, S> scalar = gram(start,start+1).inverse();
  basis.multiply_row(start+1, scalar);
  gram.multiply_col(start+1, scalar);
  gram.multiply_row(start+1, scalar);

  // Determine the appropriate scalar for clearing out the (1,1)-entry.
  if (p == 2)
    scalar = this->evaluate(basis[start+1]);
  else {
    FpElement<R,S> two(GF,2);
    scalar = -gram(start+1,start+1) / two;
  }

  // Clear the (1,1)-entry in the Gram matrix.
  basis.add_row(start+1, start, scalar);
  gram.add_col(start+1, start, scalar);
  gram.add_row(start+1, start, scalar);

  // Clear the remaining entries in the Gram matrix.
  for (size_t i = start+2; i < n; i++) {
    // Clear first row/column.
    scalar = -gram(start, i);
    basis.add_row(i, start+1, scalar);
    gram.add_col(i, start+1, scalar);
    gram.add_row(i, start+1, scalar);

    // Clear second row/column.
    scalar = -gram(start+1, i);
    basis.add_row(i, start, scalar);
    gram.add_col(i, start, scalar);
    gram.add_row(i, start, scalar);
  }

#ifdef DEBUG
  SquareMatrixFp<R, S, n> tmp = basis * original_gram;
  tmp = tmp * basis.transpose();
  assert(tmp == gram);
#endif

  // In characteristic 2, we need to recover the diagonal entries by
  //  evaluating the basis via the quadratic form.
  
  if (p == 2)
    for (size_t i = start; i < n; i++)
      gram(i, i) = this->evaluate(basis[i]);

  // basis = basis * orig_basis;
  
  return;
}

template<typename R, typename S, size_t n>
void
QuadFormFp<R, S, n>::hyperbolize_form(SquareMatrixFp<R, S, n> & gram,
				      SquareMatrixFp<R, S, n> & basis,
				      bool deterministic,
				      size_t start) const
{
  VectorFp<R, S, n> vec(this->GF);
  bool found = this->isotropic_vector(vec, start, deterministic);
  size_t dim = n - start;
  
  // The space is anisotropic.
  if (!found) {
    //    basis.set_identity();
    SquareMatrixFp<R,S,n> originalGram = gram;
    FpElement<R,S> scalar;
    if (dim == 1) {
      // Check if the (0,0)-entry is a square.
      if (gram(start,start).is_square()) {
	// If so, make it a 1.
	scalar = gram(start,start).sqrt().inverse();
	basis.multiply_row(start, scalar);
	gram.multiply_col(start, scalar);
	gram.multiply_row(start, scalar);
      }
      return;
    }
    if (this->GF->prime() == 2) {
      // Make the (0,0)-entry equal to 1.
#ifdef DEBUG
      assert(gram(start,start).is_square());
#endif
      scalar = gram(start,start).sqrt().inverse();
      basis.multiply_row(start, scalar);
      gram.multiply_col(start, scalar);
      gram.multiply_row(start, scalar);

      // Make the (0,1)-entry equal to 1.
      scalar = gram(start,start+1).inverse();
      basis.multiply_row(start+1, scalar);
      gram.multiply_col(start+1, scalar);
      gram.multiply_row(start+1, scalar);

      return;
    }
    
    // Clear the (0,1)-entry.
    scalar = -gram(start,start+1)/gram(start,start);
    basis.add_row(start+1, start, scalar);
    gram.add_col(start+1, start, scalar);
    gram.add_row(start+1, start, scalar);

    // If the (1,1)-entry is a square, make it the first entry.
    if (gram(start+1,start+1).is_square()) {
      basis.swap_rows(start,start+1);
      gram.swap_rows(start,start+1);
      gram.swap_cols(start,start+1);
    }

    bool is_square[2];
    for (size_t i = 0; i < 2; i++) {
      // Check if the (i,i)-entry is a square then clear it, if so.
      is_square[i] = gram(start+i,start+i).is_square();
      if (is_square[i]) {
	scalar = gram(start+i,start+i).sqrt().inverse();
	basis.multiply_row(start+i, scalar);
	gram.multiply_col(start+i, scalar);
	gram.multiply_row(start+i, scalar);
      }
    }

    // If neither are squares, make them -1 (note that this occurs
    //  if and only if -1 is not a square).
    if ((!is_square[0]) && (!is_square[1])) {
      for (size_t i = 0; i < 2; i++) {
#ifdef DEBUG
	assert((-gram(start+i,start+i)).is_square());
#endif
	scalar = (-gram(start+i,start+i)).sqrt().inverse();
	basis.multiply_row(start+i, scalar);
	gram.multiply_col(start+i, scalar);
	gram.multiply_row(start+i, scalar);
      }
    }

    return;
  }

#ifdef DEBUG
  assert(this->evaluate(vec) == 0);
#endif

  // Attempt to split a hyperbolic plane from the form.
  split_hyperbolic_plane(vec, gram, basis, start);

  // Determine how many dimensions we need to split off.
  size_t lower_dim = gram[0].is_zero() ? 1 : 2;

  if (dim > lower_dim) {
    // Split the hyperbolic plane from the form.
    QuadFormFp<R, S, n> q_split(gram);
    // !! TODO - check maybe we have to replace basis here
    SquareMatrixFp<R, S, n> newbasis(this->GF);
    newbasis.set_identity();
    q_split.hyperbolize_form(gram, newbasis, deterministic, start + lower_dim);
    basis = newbasis * basis;
  }

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "After hyperbolize_form with start = " << start << "." << std::endl;
  std::cerr << "Resulting gram matrix is " << std::endl;
  gram.pretty_print(std::cerr);
  std::cerr << ", ";
  std::cerr << "Resulting basis is " << std::endl;
  basis.pretty_print(std::cerr);
  std::cerr << std::endl;
#endif 
  
  return;
}

template<typename R, typename S, size_t n>
void QuadFormFp<R, S, n>::decompose(SquareMatrixFp<R, S, n> & gram,
				    SquareMatrixFp<R, S, n> & basis,
				    bool deterministic) const
{
  basis.set_identity();
  hyperbolize_form(gram, basis, deterministic);

#ifdef DEBUG_LEVEL_FULL
  std::cerr << "After hyperbolize_form." << std::endl;
  std::cerr << "Resulting gram matrix is " << gram << ", ";
  std::cerr << "Resulting basis is " << basis << std::endl;
  
  // Verify that everyhing we've done is correct.
  SquareMatrixFp<R, S, n> temp1(this->GF);
  SquareMatrixFp<R, S, n> temp2(this->GF);
  if (this->GF->prime() == 2) {
    // Verify that the basis evaluates correctly on the form.
    for (size_t i = 0; i < n; i++) {
      assert(this->evaluate(basis[i]) == gram(i,i));
    }
    // Zero out the diagonal to verify the bilinear form is correct.
    // Recall that in characteristic 2, the associated bilinear form
    // is b_q(x,y) = q(x+y) - q(x) - q(y), with zeros on the diagonal
    temp1 = gram;
    temp2 = this->bilinear_form();
    for (size_t i = 0; i < n; i++) {
      temp1(i,i) = 0;
      temp2(i,i) = 0;
    }
  }
  else {
    temp1 = gram;
    temp2 = this->bilinear_form();
  }
  // Verify that the bilinear forms are similar.
  assert(basis * temp2 * basis.transpose() == temp1);
#endif

  // Let's bubble the basis vectors which belong to the radical to the
  //  end of the basis list.
  size_t rad = 0;
  size_t pos = n;
  while (pos >= 1) {
    if (gram[pos-1].is_zero()) {
      rad++;
      for (size_t i = pos; i < n; i++) {
	basis.swap_rows(i-1,i);
	gram.swap_rows(i-1,i);
	gram.swap_cols(i-1,i);
      }
    }
    pos--;
  }
  // Let's put the hyperbolic planes in our standard antidiagonal form.

  // The upper index of the hyperbolic space.
  size_t upper = n + 1 - rad;
  do {
    upper--;
  } while ((upper >= 1) && (gram(upper-1,upper-1) != 0));

  // Indices of the basis vectors we'll be swapping.
  size_t i = 1;
  size_t j = upper;
  // Keep swapping basis vectors until j is less than or equal to i. Note
  //  that if there are no hyperbolic planes, this does nothing.
  while (i < j) {
    basis.swap_rows(i-1,j-2);
    gram.swap_cols(i-1,j-2);
    gram.swap_rows(i-1,j-2);
    i += 2;
    j -= 2;
  }
  // Since we did everything with row vectors, we need to transpose the
  //  basis, so that the subsequent code that utilizes it doesn't break.
  basis = basis.transpose();
  return;
}

// !! TODO - all F2 operations (including this one) can be made faster
template<typename R, typename S, size_t n>
FpElement<R, S> QuadFormFp<R, S, n>::evaluate_p2(const VectorFp<R, S, n>& v)
  const
{
  FpElement<R, S> val(this->GF, 0);
  for (size_t i = 0; i < n; i++) {
    val += this->B_Fp(i, i) * v[i];
    for (size_t j = i+1; j < n; j++)
      val += this->B_Fp(i,j) * v[i] * v[j];
  }
  return val;
}

template<typename R, size_t n>
std::vector<uint8_t>
QuadForm_Base<R,n>::bit_transpose(const std::vector< uint8_t > & mat)
{
#ifdef DEBUG
  assert(n+1 <= 8);
#endif
  std::vector<uint8_t> trans(n+1);

  for (uint8_t row = 0; row <= n; row++) {
    trans[row] = 0;
    for (uint8_t col = 0; col < mat.size(); col++) {
      trans[row] |= (((mat[col] >> row) & 1) << col);
    }
  }

  return trans;
}

// returns the transformation and the rank
// performs the echelonization in place
template<typename R, size_t n>
uint8_t QuadForm_Base<R,n>::bit_echelon_form(std::vector< uint8_t > & mat,
					     std::vector< uint8_t > & trans)
{
#ifdef DEBUG
  assert(mat.size() <= 8);
#endif
  
  trans.resize(mat.size());
  
  for (size_t row = 0; row < mat.size(); row++)
    trans[row] = (1 << row);
  
  uint8_t pivot_row;
  pivot_row = 0;
  uint8_t pivot_col;
  pivot_col = 0;
 
  uint8_t row;  
  uint8_t val;
  
  while ((pivot_row < mat.size()) && (pivot_col <= n)) {
    val = 0;
    for (row = pivot_row ; (!val) && (row < mat.size()); row++) {
      val = (mat[row] >> pivot_col) & 1;
    }
    if (!val) {
      pivot_col++;
    }
    else {
      row--;
      if (row != pivot_row) {
	// swapping rows
	mat[pivot_row] ^= mat[row];
	mat[row] ^= mat[pivot_row];
	mat[pivot_row] ^= mat[row];
      
	trans[pivot_row] ^= trans[row];
	trans[row] ^= trans[pivot_row];
	trans[pivot_row] ^= trans[row];
      }
      
      for (row = pivot_row+1; row < mat.size(); row++) {
	val = (mat[row] >> pivot_col) & 1;
	if (val) {
	  mat[row] ^= mat[pivot_row];
	  trans[row] ^= trans[pivot_row];
	}
      }
      
      pivot_row++;
      pivot_col++;
    }
  }
  return pivot_row;
}

template<typename R, size_t n>
std::vector<uint8_t>
QuadForm_Base<R,n>::kernel(const std::vector< uint8_t > & mat)
{
  std::vector<uint8_t> ker;

  std::vector<uint8_t> mat_t = bit_transpose(mat);
  
  uint8_t rank = bit_echelon_form(mat_t, ker);
  // getting the zero rows

  ker.erase(ker.begin(), ker.begin() + rank);
  
  return ker;
}
  
// !! TODO - use bit slicing to make this faster
// Also does not need to compute the rank every time

template<typename R, size_t n>
bool QuadForm_Base<R,n>::sign_normalization_fast(SquareMatrix<R, n> & qf,
						 Isometry<R,n> & isom)
{
  bool is_reduced = true;

  // This assumes n < 8
#ifdef DEBUG
  assert(n < 8);
#endif
  std::vector< uint8_t > bb_vecs;
  uint8_t vec;
  // vec after reduction - for echelon form
  uint8_t ech_vec;
  
  // There should be a more efficient way of doing this,
  // but it will help me keep track of things for now
  // save the pivots of each row, this is always sorted
  std::vector<uint8_t> pivots;
  
  // for each k (from 0 to n-1) save the row number
  // in which we will want it to be placed
  // If it is already a pivot, this will be -1
  int8_t place_pivots[n] = {0};

  // the position of the row where k is s pivot.
  // If it is not, then it is -1.
  int8_t inv_pivots[n];
  for (uint8_t i = 0; i < n-1; i++) inv_pivots[i] = -1;
  
  for (size_t j = 1; j < n; j++) {
    // vec will always have the only the bits k and k+j on 
    vec = 1 | (1 << j);
    for (size_t k = 0; k < n-j; k++) {
      if (qf(k,k+j) != 0) {
	int lead = k;
	ech_vec = vec;

	if (qf(k,k+j) < 0)
	  ech_vec |= (1 << n);

	// while we already have this as pivot, we 
	while ((lead >= 0) && (inv_pivots[lead] >= 0) && (lead != n)) {
	  ech_vec ^= bb_vecs[inv_pivots[lead]];
	  lead = ffs(ech_vec)-1;
	}

	// If it is not a pivot, we put it in its proper place
	// and update the arrays tracking the pivots
	if ((lead >= 0) && (lead != n)) {
	  // echelonize the rows already in bb_vecs above the pivot
	  for (uint8_t row = 0; row < place_pivots[lead]; row++) {
	    uint8_t bit = (bb_vecs[row] >> lead) & 1;
	    if (bit) bb_vecs[row]^= ech_vec;
	  }
	  bb_vecs.insert(bb_vecs.begin()+place_pivots[lead], ech_vec);
	  inv_pivots[lead] = place_pivots[lead];
	  place_pivots[lead] = -1;
	  // next time we will get any of these, we will want
	  // them after this one
	  for (size_t r = lead+1; r < n; r++) {
	    place_pivots[r]++;
	    if (inv_pivots[r] >= 0)
	      inv_pivots[r]++;
	  } 
	}
      }
      vec <<= 1;
    }
  }

  // The last row of ker should now be a solution to the affine equation
  // The rows above are the kernel
  std::vector<uint8_t> ker_bit = kernel(bb_vecs);

#ifdef DEBUG_LEVEL_FULL
  Isometry<R, n> s;
  for (size_t i = 0; i < n; i++)
    s(i,i) = ((ker_bit[ker_bit.size()-1] >> i) & 1) ? -1 : 1;
  if (s.transform(qf) == qf) {
    is_reduced = true;
    // to be compatible with magma implementation for debugging
    for (size_t i = 0; i < n; i++) s(i,i) = 1;
  }
  SquareMatrix<R, n> qf_orig = qf;
  Isometry<R,n> isom_orig = isom;
#endif 
  is_reduced = true;
  
  if (!ker_bit.empty()) {
    for (size_t row = 0; row+1 < n; row++) {
      uint8_t bit_row = (ker_bit[ker_bit.size()-1] >> row) & 1;
      for (size_t col = row+1; col < n; col++) {
	uint8_t bit_col = (ker_bit[ker_bit.size()-1] >> col) & 1;
	if (bit_col^bit_row) {
	  qf(row,col) = -qf(row,col);
	  qf(col,row) = -qf(col,row);
	  is_reduced = false;
	}
      }
    }
#ifdef DEBUG_LEVEL_FULL
    assert(s.transform(qf_orig) == qf);
#endif
    if (!is_reduced)
      for (size_t col = 0; col < n; col++)
	if ((ker_bit[ker_bit.size()-1] >> col) & 1) {
	  for (size_t row = 0; row < n; row++)
	    isom(row,col) = -isom(row,col);
	}
#ifdef DEBUG_LEVEL_FULL
    assert(isom_orig*s == isom);
#endif
  }
  
  return is_reduced;
}


template<typename R, size_t n>
bool
QuadForm_Base<R,n>::sign_normalization(SquareMatrix<R, n> & qf,
				       Isometry<R,n> & isom,
				       std::set< Isometry<R, n> > & auts,
				       bool calc_aut)
{
  if (!calc_aut)
    return sign_normalization_fast(qf, isom);
  return sign_normalization_slow(qf, isom, auts);
}
