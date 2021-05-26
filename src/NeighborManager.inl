#include "Polynomial.h"

template<typename R, typename S, typename T, size_t n>
NeighborManager<R,S,T,n>::NeighborManager(const QuadForm<T, n>& q,
					  std::shared_ptr<Fp<R,S>> GF,
					  size_t k)
  : vec(GF)
{
  T p = GF->prime();
  
  this->q = q;
  this->disc = q.discriminant();

  std::shared_ptr<QuadFormFp<R,S,n> > qp = q.mod(GF);

  this->b = qp->bilinear_form();
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      this->quot_gram(i,j) = this->q.bilinear_form()(i,j) % (p*p);
	
  this->GF = GF;
  assert(qp->isotropic_vector(this->vec));

#ifdef DEBUG
  R prime = GF->prime();
  if (prime != 2) assert( qp->evaluate(vec) == 0 );
#endif

  this->p_std_gram = std::make_shared<SquareMatrixFp<R, S, n> >(GF);
  this->p_basis = std::make_shared<SquareMatrixFp<R, S, n> >(GF);

#ifdef DEBUG
  qp->decompose(*p_std_gram, *p_basis, true);
#else
  qp->decompose(*p_std_gram, *p_basis);
#endif
  
  this->p_q_std = std::make_shared<PolynomialFp<R,S> >(*p_std_gram);

#ifdef DEBUG
  std::cerr << "Performed Witt Decomposition on" << std::endl;
  qp->bilinear_form().pretty_print(std::cerr);
  std::cerr << "Resulting gram matrix is " << std::endl;
  p_std_gram->pretty_print(std::cerr);
  std::cerr << "Resulting basis is " << std::endl;
  p_basis->pretty_print(std::cerr);
#endif

  // Count the rows at the end of the matrix which are exactly zero.
  size_t idx = n;
  while ((idx >= 1) && (*p_std_gram)[idx-1].is_zero()) {
    idx--;
  }

  // The dimension of the radical.
  this->rad_dim = n - idx;

  // Determine the dimension of the totally hyperbolic subspace.
  idx = 1;
  while ((idx <= n-rad_dim) && ((*p_std_gram)(idx-1, idx-1) == 0) )
    idx++;

  // Dimension of the anistotropic subspace.
  this->aniso_dim = n - rad_dim - idx + 1;

  // The number of hyperbolic planes in the Witt decomposition.
  this->witt_index = (idx - 1) / 2;

  this->pivots = __pivots(n-rad_dim, aniso_dim, k);
  this->pivot_ptr = 0;
  this->k = k;
  this->skew_dim = k*(k-1)/2;
  this->p_skew = std::make_shared< MatrixFp<R,S> >(this->GF, k, k);
}

//!! TODO - make gram work only modulo p^2

template<typename R, typename S, typename T, size_t n>
SquareMatrix<T,n>
NeighborManager<R,S,T,n>::__gram(const SquareMatrix<T, n> & B, bool quot) const
{
  T p = this->GF->prime();
  SquareMatrix<T, n> gram;
  
  if (p == 2)
    gram = B * this->q.bilinear_form() * B.transpose();
  else
    gram =  B * (this->quot_gram) * B.transpose();
  
  // !! TODO - this isn't necessary only for debugging versus magma
  if ((quot) || (p != 2))
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++)
	gram(i,j) = gram(i,j) % (p*p);
  
  return gram;
}

template<typename R, typename S, typename T, size_t n>
void NeighborManager<R,S,T,n>::lift_subspace()
{
  T p = this->GF->prime();

#ifdef DEBUG
  assert(this->pivot_ptr >= 1);
#endif
  
  // Get the pivots for the bases of the isotropic subspaces.
  std::vector<size_t> pivots = this->pivots[this->pivot_ptr-1];

#ifdef DEBUG
  std::cerr << "before lifting, p_basis is " << std::endl << (*this->p_basis);
  std::cerr << std::endl;
  std::cerr << "iso_subspace is " << this->iso_subspace << std::endl;
  std::cerr << "pivots = ";
  for (size_t i = 0; i < pivots.size(); i++)
    std::cerr << pivots[i];
  std::cerr << std::endl;
#endif
  
  // Set up the correct basis vectors.
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = pivots[i]+1; j < n; j++)
      this->p_basis->add_col(pivots[i], j, this->iso_subspace[i][j]);

#ifdef DEBUG
  std::cerr << "the correct basis vectors are" << std::endl << (*this->p_basis);
  std::cerr << std::endl;
#endif
  
  // Extract our target isotropic subspace modulo p
  std::vector< VectorFp<R, S, n> > x,z,u;
  for (size_t i = 0; i < this->k; i++) {
    x.push_back(p_basis->transpose()[pivots[i]]);
  }

#ifdef DEBUG
  std::cerr << "x = " << x << std::endl;
#endif
  
  // Extract the hyperbolic complement modulo pR.
  std::vector<size_t> paired(this->k);
  size_t h_dim = 2 * this->witt_index; 
  for (size_t i = 0; i < this->k; i++)
    paired[i] = h_dim - 1 - pivots[this->k-1-i];
  for (size_t i = 0; i < this->k; i++) {
    z.push_back(p_basis->transpose()[paired[i]]);
  }

#ifdef DEBUG
  std::cerr << "z = " << z << std::endl;
#endif
  
  // Extract the remaining basis vectors.
  std::set<size_t> exclude;
  exclude.insert(pivots.begin(), pivots.end());
  exclude.insert(paired.begin(), paired.end());
  for (size_t i = 0; i < n; i++) {
    std::set<size_t>::const_iterator iter = exclude.find(i);
    if (iter == exclude.end())
      u.push_back(p_basis->transpose()[i]);
  }

#ifdef DEBUG
  std::cerr << "u = " << u << std::endl;
#endif
  
  // Convert to coordinates modulo p^2.
  X.resize(this->k);
  Z.resize(this->k);
  U.resize(n - 2*this->k);
  
  // Build the coordinate matrix.
  // !! TODO - the mod p is not necessary, good for debugging
  SquareMatrix<T, n> B;
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(i,j) = X[i][j] = x[i][j].lift() % p;

  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(this->k+i,j) = Z[i][j] = z[i][j].lift() % p;

  for (size_t i = 0; i < n - 2*this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(2*this->k+i,j) = U[i][j] = u[i][j].lift() % p;

#ifdef DEBUG
  std::cerr << "X = " << X << std::endl;
  std::cerr << "Z = " << Z << std::endl;
  std::cerr << "U = " << U << std::endl;
#endif
  
  // Compute the Gram matrix of the subspace with respect to the spaces
  //  we will perform the following computations upon.

  SquareMatrix<T,n> gram = __gram(B);

  // Lift Z so that it is in a hyperbolic pair with X modulo p^2.
  std::vector< Vector<T, n> > Z_new(k);
  for (size_t i = 0; i < this->k; i++) {
    Z_new[i] = Z[i];
    for (size_t j = 0; j < this->k; j++) {
      T delta = (i == j) ? 1 : 0;
      // we split the operation due to signed type issues
      T a = gram(this->k-j-1, i + this->k);
      // a nonnegative value with the same residue mod p*p
      // !!! TODO - we might be able to get rid of that
      // since T is always signed - check!
      a = (a / (p*p) + 1)*p*p-a+delta;
      if (a >= (p*p))
	a -= p*p;
      Z_new[i] += a * Z[j];
    }
  }
  Z = Z_new;
  
#ifdef DEBUG
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      Z[i][j] = Z[i][j] % (p*p);
  
  std::cerr << "after setting <X,Z> = 1" << std::endl;
  std::cerr << "Z = " << Z << std::endl;

  // Verify that X and Z form a hyperbolic pair.
  // Compute the Gram matrix thusfar.
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(this->k+i,j) = Z[i][j];
  
  SquareMatrix<T, n> temp = __gram(B);
  for (size_t i = 0; i < k; i++)
    for (size_t j = 0; j < k; j++)
      assert(temp(i, k+j) % (p*p) == ((i+j == k-1) ? 1 : 0));	
#endif
  
  if (p == 2) {
    for (size_t i = 0; i < this->k; i++)
      for (size_t j = 0; j < n; j++)
	B(this->k+i,j) = Z[i][j];
    gram = __gram(B, false);
  }
  // Lift X so that it is isotropic modulo p^2.
  std::vector< Vector<T, n> > X_new(k);
  T half = (p*p+1)/2;
  for (size_t i = 0; i < this->k; i++) {
    X_new[i] = X[i];
    T gram2 = gram(i,i)/2 + ((gram(i,i) % 2 == 0) ? 0 : half);
    for (size_t j = this->k-1-i; j < this->k; j++) {
      T scalar = (i+j == k-1) ? gram2 : gram(i, this->k-1-j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      X_new[i] +=  scalar  * Z[j];
    }
  }
  X = X_new;

#ifdef DEBUG
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      X[i][j] = X[i][j] % (p*p);
  std::cerr << "after setting <X,X> = 0" << std::endl;
  std::cerr << "X = " << X << std::endl;
  
  // Verify that X is isotropic modulo p^2.
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(i,j) = X[i][j];

  // The Gram matrix on this basis.
  temp = __gram(B);

  // Verify all is well.
  for (size_t i = 0; i < k; i++)
    for (size_t j = 0; j < k; j++)
      assert(temp(i,j) % (p*p) == 0);
  
#endif

  // Lift Z so that it is isotropic modulo p^2.
  for (size_t i = 0; i < this->k; i++) {
    Z_new[i] = Z[i];
    for (size_t j = this->k-1-i; j < this->k; j++) {
      T scalar = (i+j == k-1) ? 2 : 1;
      scalar = gram(this->k+i, 2*this->k-1-j) / scalar;
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      Z_new[i] += scalar * X[j];
    }
  }
  Z = Z_new;

#ifdef DEBUG
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      Z[i][j] = Z[i][j] % (p*p);
  std::cerr << "after setting <Z,Z> = 0" << std::endl;
  std::cerr << "Z = " << Z << std::endl;
  
  // Verify that Z is isotropic modulo p^2.
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(this->k+i,j) = Z[i][j];

  // The Gram matrix on this basis.
  temp = __gram(B);

  // Verify all is well.
  for (size_t i = 0; i < k; i++)
    for (size_t j = 0; j < k; j++)
      assert(temp(k+i,k+j) % (p*p) == 0);
  
#endif
  
  // The Gram matrix thusfar.
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(i,j) = X[i][j];

  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(this->k+i,j) = Z[i][j];

  gram = __gram(B);

  // Make U orthogonal to X+Z.
  for (size_t i = 0; i < k; i++)
    for (size_t j = 0; j < n - 2*k; j++) {
      // Clear components corresponding to X.
      T scalar = gram(2*k-1-i, 2*k+j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      U[j] += scalar * X[i];
      
      // Clear components corresponding to Z.
      scalar = gram(k-1-i, 2*k+j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      U[j] += scalar * Z[i];
    }

#ifdef DEBUG
  for (size_t i = 0; i < n-2*this->k; i++)
    for (size_t j = 0; j < n; j++)
      U[i][j] = U[i][j] % (p*p);
  std::cerr << "after setting <U,X+Z> = 0" << std::endl;
  std::cerr << "U = " << U << std::endl;
  
  // Verify that U is now orthogonal to X+Z.
  for (size_t i = 0; i < n-2*this->k; i++)
    for (size_t j = 0; j < n; j++)
      B(2*this->k+i,j) = U[i][j];

  // The Gram matrix on this basis.
  temp = __gram(B);

  // Verify all is well.
  for (size_t i = 0; i < 2*k; i++)
    for (size_t j = 2*k; j < n; j++)
      assert(temp(i,j) % (p*p) == 0);

  // make sure that all the entries of U are between 0 and p^2
  for (size_t i = 0; i < n-2*k; i++)
    for (size_t j = 0; j < n; j++)
      U[i][j] = U[i][j] % (p*p);
#endif

  return;
}

template<typename R, typename S, typename T, size_t n>
void NeighborManager<R,S,T,n>::next_isotropic_subspace()
{
  if (this->params.empty()) {
    // Move to the next pivot.
    this->pivot_ptr++;
    
    // If we've exceeded the list of pivots, we're done.
    if (this->pivot_ptr > this->pivots.size()) {
      // Reset the pivot pointer so that we can loop through
      //  the isotropic subspaces again if needed.
      this->pivot_ptr = 0;
      this->iso_subspace.clear();
      return;
    }
    
    // Initialize the new pivot.
    __initialize_pivots();
  }

  // The list of evaluation values.
  FpElement<R,S> zero(this->GF, 0);
  std::vector<FpElement<R,S> > eval_list(n*k, zero);

  // Produce the isotropic subspace corresponding to the current
  //  parameters.
  for (size_t i = 0; i < this->params.size(); i++)
    eval_list[this->free_vars[i]] = this->params[i];

  // The basis for the current isotropic subspace.
  for (size_t i = 0; i < this->k; i++) {
    VectorFp<R,S,n> vec_fp(this->GF);
    for (size_t j = 0; j < n; j++) {
      vec_fp[j] = (*(this->p_isotropic_param))(i,j).evaluate(eval_list);
    }
    this->iso_subspace.push_back(vec_fp);
  }

  if (this->free_vars.size() != 0) {
    // The current position in the parameterization.
    size_t pos = 0;
    // Terminate loop once we found the next new subspace, or we
    //  hit the end of the list.
    do {
      // Increment position.
      pos++;
      // Manually move to the next element.
      this->params[pos]++;
    } while ((pos != this->free_vars.size()) && (this->params[pos] == 0));
  }

  // If we've hit the end of the list, indicate we need to move on to the
  //  next pivot.

  bool all_zero = true;
  for (size_t i = 0; i < this->params.size(); i++)
    if (this->params[i] != 0) {
      all_zero = false;
      break;
    }

  if (all_zero) {
    this->params.clear();
  }
  
  return;
}

template<typename R, typename S, typename T, size_t n>
inline GenusRep<T, n> NeighborManager<R,S,T,n>::get_reduced_neighbor_rep(R t) const
{
  GenusRep<T, n> rep;
  rep.q = this->get_neighbor(t, rep.s);
  rep.q = QuadForm<T, n>::reduce(rep.q, rep.s);
  return rep;
}

// to representative of the line
template<typename R, typename S, typename T, size_t n>
Vector<R, n> NeighborManager<R,S,T,n>::transform_vector(const GenusRep<T, n>& dst,
							Vector<R, n> src)
{
  Vector<T, n> temp;
  VectorFp<R, S, n> temp_mod = GF->mod(src);
  for (size_t i = 0; i < n; i++)
    temp[i] = temp_mod[i].lift();

  R p = GF->prime();

  // should that be inverse mod p?
  Isometry<T, n> sinv = dst.s.inverse();
  temp = sinv * temp;

#ifdef DEBUG
  for (size_t i = 0; i < n; i++)
    assert( temp[i] % p == 0 );
#endif

  for (size_t i = 0; i < n; i++)
    temp[i] /= p;

  Vector<R, n> vec;
  for (size_t i = 0; i < n; i++)
    vec[i] = GF->mod(temp[i]).lift();

  for (size_t i = n; i > 0; i--) {
    if (vec[i-1] != 0)
      {
	R inv = GF->inverse(vec[i-1]);
	for (size_t j = 0; j < i-1; j++)
	  vec[j] = GF->mod(GF->mul(vec[j], inv)).lift();
	vec[i-1] = 1;
	break;
      }
  }

  return vec;
}

template<typename R, typename S, typename T, size_t n>
void NeighborManager<R,S,T,n>::update_skew_matrix(size_t & row, size_t & col)
{
  bool done;
  do {
    // Flag for determining whether we are done updating
    //  the skew matrix.
    done = true;
     
    // Increment value of the (row,col) position.
    (*(this->p_skew))(row, col)++;
      
    // Update the coefficient of the skew matrix reflected
    //  across the anti-diagonal.
    (*(this->p_skew))(k-1-col, k-1-row) = -(*(this->p_skew))(row,col);
      
    // If we've rolled over, move on to the next position.
    if ((*(this->p_skew))(row,col) == 0) {
      // The next column of our skew matrix.
      col++;
      // Are we at the end of the column?
      if (row+col == k-1) {
	// Yes. Move to the next row.
	row++;
	// And reset the column.
	col = 0;
      }
      // Indicate we should repeat another iteration.
      done = false;
    }
  } while ((!done) && (row+col != k-1));
  return;
}

template<typename R, typename S, typename T, size_t n>
void NeighborManager<R,S,T,n>::update_skew_space()
{
  T p = this->GF->prime();
  // Update the skew space.
  for (size_t i = 0; i < k ; i++) {
    for (size_t j = 0; j < k; j++){
      // !! TODO - I got rid here of X_skew,
      // check that it sisn't destroy anything
      T val = (*(this->p_skew))(i,j).lift();
      this->X[i] += p * (val * this->Z[j]);
    }
  }
}

template<typename R, typename S, typename T, size_t n>
void NeighborManager<R,S,T,n>::get_next_neighbor(void)
{
  size_t row,col;
  // The starting position of the skew vector to update.
  row = 0;
  col = 0;
  // Update the skew matrix (only for k >= 2).
  if (this->skew_dim != 0) {
    this->update_skew_matrix(row, col);
  }
  
  // If we haven't rolled over, update the skew space and return...
  if (row+col < k-1) {
    this->update_skew_space();
    return;
  }

  // ...otherwise, get the next isotropic subspace modulo p.
  this->next_isotropic_subspace();

  // Lift the subspace if we haven't reached the end of the list.
  if (!(this->iso_subspace.empty())) {
    this->lift_subspace();
  }
  return;
}

template<typename R, typename S, typename T, size_t n>
QuadForm<T, n> NeighborManager<R,S,T,n>::build_neighbor(Isometry<T, n>& s) const
{
  T p = GF->prime();
  T p2 = p*p;
  T p3 = p2*p;
  SquareMatrix<T, n> qq;

  // fill the isomtery by the X,Z,U
  // if lift_subspace was successful,
  // <X,X>, <X,U>,<Z,Z>,<Z,U> in p^2 and <X,Z> = 1 mod p^2

  // For now, we follow thw magma implementation
  // start by scaling the basis
  
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      s(i,j) = this->X[i][j];

  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < n; j++)
      s(this->k+i,j) = p2*this->Z[i][j];

  for (size_t i = 0; i < n-2*this->k; i++)
    for (size_t j = 0; j < n; j++)
      s(2*this->k+i,j) = p*this->U[i][j];

  SquareMatrix<T,n> hermite = s.a.hermite_form(p3);
  
  s = hermite.transpose();
	 
  //	s.swap_cols(0, pivot);
  // need to adjust determinant for s to be in SO
  // This transforms using the isometry
  qq = s.transform(q.bilinear_form(), 1);

#ifdef DEBUG
  assert( qq(0,0) % p == 0 );
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      assert(qq(i,j) % p2 == 0);
#endif

  // we have to rescale
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      qq(i,j) /= p2;
  
  /*
  // Here we simulate (X,Z,U) |-> (X/p,pZ,U)
  // since we don't want to introduce denominators
  
  // 1. divide <X,X> by p^2
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < this->k; j++)
      qq(i,j) /= (p*p);

  // 2. multiply <Z,Z> by p^2
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < this->k; j++)
      qq(k+i,k+j) *= (p*p);

  // 3. divide <X,U> by p
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 2*this->k; j < n; j++) {
      qq(i,j) /= p;
      qq(j,i) /= p;
    }

  // 4. multiply <Z,U> by p
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 2*this->k; j < n; j++) {
      qq(this->k+i,j) *= p;
      qq(j, this->k+i) *= p;
    }
  */
  QuadForm<T, n> retval(qq);
	
  if (std::is_same<T,Z64>::value)
    {
      // If we're not using arbitrary precision, throw an exception if
      // the discriminant of the p-neighbor isn't correct.
      if (retval.discriminant() != this->disc)
	{
	  throw std::overflow_error(
				    "An overflow has occurred. The p-neighbor's discriminant "
				    "does not match the original.");
	}
    }
  return retval;
}

template<typename R, typename S, typename T, size_t n>
Vector<R, n> NeighborManager<R,S,T,n>::isotropic_vector_p2(R t) const
{
  Vector<R, n> res;
  res[n-1] = 1;
    
  // Stub
  // !! TODO - do something appropriate here
       
  return res;
}

// A helper function for computing valid pivots.
template<typename R, typename S, typename T, size_t n>
std::vector< std::vector<size_t> >
NeighborManager<R,S,T,n>::__pivots(size_t dim, size_t aniso, size_t k)
{
  std::vector< std::vector<size_t> > pivs;
  // Base case.
  if (k == 1) {
    for (size_t i = 0; i < dim - aniso; i++) {
      std::vector<size_t> singleton(1);
      singleton[0] = i;
      pivs.push_back(singleton);
    }
    return pivs;
  }

  // Retrieve lower-dimensional maximal pivots.
  pivs = __pivots(dim-2, aniso, k-1);
  for (size_t i = 0; i < pivs.size(); i++)
    for (size_t j = 0; j < pivs[i].size(); j++)
      pivs[i][j]++;

  size_t num = pivs.size();
  // Determine the first set of pivots.
  pivs.insert(pivs.end(), pivs.begin(), pivs.end());
  for (size_t i = 0; i < num; i++){
    pivs[i].insert(pivs[i].begin(), 0);
    pivs[i+num].push_back(dim-aniso-1);
  } 

  // Add additional pivots when we're not in the maximal case.
  if (2*k <= dim - aniso) {
    std::vector< std::vector<size_t> > pivs2 = __pivots(dim-2, aniso, k);
    for (size_t i = 0; i < pivs2.size(); i++)
      for (size_t j = 0; j < pivs2[i].size(); j++)
	pivs2[i][j]++;
    pivs.insert(pivs.begin(), pivs2.begin(), pivs2.end());
  }
  return pivs;
}

template<typename R, typename S, typename T, size_t n>
void NeighborManager<R,S,T,n>::__initialize_pivots(void)
{
#ifdef DEBUG
  assert(this->pivot_ptr >= 1);
#endif 
  std::vector<size_t> pivot = this->pivots[this->pivot_ptr-1];
  size_t rank = (this->k)*n;
  
  // Keep a list of non-free variables from which we will determine the
  //  free variables when we are done.
  std::vector<size_t> remove;

  // Initialize matrix that will determine parameterization.
  std::vector< PolynomialFp<R, S> > data;
  for (size_t i = 0; i < rank; i++) {
    PolynomialFp<R,S> x_i(this->GF, i);
    data.push_back(x_i);
  }
  this->p_isotropic_param =
    std::make_shared< Matrix< PolynomialFp<R, S> > >(data, this->k, n);

  FpElement<R,S> zero(this->GF, 0);
  FpElement<R,S> one(this->GF, 1);
  // Setup the columns corresponding to the pivots.
  for (size_t row = 0; row < this->k; row++)
    for (size_t col = 0; col < this->k; col++) {
      (*p_isotropic_param)(row, pivot[col]) = (row == col) ? one : zero;
      remove.push_back(row*n + pivot[col]);
    }

  // Clear the rows prior to the pivot positions (but not the radical).
  for (size_t row = 0; row < this->k; row++)
    for (size_t col = 0; col < pivot[row]; col++) {
      (*p_isotropic_param)(row, col) = zero;
      remove.push_back(row*n + col);
    }

  // Check if one or more of the anisotropic coordinates need to be zero.
  for (size_t row = 0; row < k; row++) {
    if (pivot[row] > this->witt_index) {
      for (size_t col = 0; col < this->aniso_dim; col++) {
	(*p_isotropic_param)(row, n-1-rad_dim-col) = zero;
	remove.push_back((row+1)*n-1-rad_dim-col);
      }
    }
  }

  // Determine the number of rows of the matrix that we'll echelonize.
  size_t rows = k*(k+1)/2;

  // Here we will save the quadratic polynomials
  data.clear();
  
  
  // The matrix that we're going to echelonize.
  MatrixFp<R, S> mat(this->GF, rows, rank);
  
  // The current row to fill in in the matrix.
  size_t row = 0;
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = i; j < this->k; j++) {
      // The appropriate vector that we want to be isotropic.
      std::vector<PolynomialFp<R, S> > vec;
      for (size_t r = 0; r < n; r++)
	vec.push_back((i == j) ? (*p_isotropic_param)(i,r) :
		      (*p_isotropic_param)(i,r) + (*p_isotropic_param)(j,r));
   
      PolynomialFp<R, S> f = p_q_std->evaluate(vec);
      // Degree 2 terms are inhomogenous.
      data.push_back(-f.quadratic_part());
      // The other terms are linear
      // so fill in the matrix accordingly.
      std::vector< FpElement<R,S> > l = f.linear_part(rank);
      for (size_t i = 0; i < rank; i++)
	mat(row, i) = l[i];
      // Move on to the next row.
      row++;
    }

  // Compute the Echelon form of the matrix.
  MatrixFp<R,S> trans(this->GF, rows, rows);
  MatrixFp<R,S>::row_echelon(mat, trans);

  // The evaluation list for replacing variables with their dependence
  //  relations.
  std::vector< PolynomialFp<R, S> > eval_list;
  for (size_t i = 0; i < rank; i++) {
    PolynomialFp<R,S> x_i(this->GF, i);
    eval_list.push_back(x_i);
  }

  for (size_t i = 0; i < rows; i++) {
    // Find the pivot in the i-th row.
    size_t c = 0;
    while ((c < rank) && (mat(i,c) != 1)) c++;
    // Add this pivot to the list of non-free variables.
    remove.push_back(c);

    // If the row is entirely zero, skip it.
    if (c >= rank) continue;

    // Build the multinomial for which x_c is dependent.
    //    PolynomialFp<R,S> f = mat(i, rank);
    PolynomialFp<R,S> f(this->GF);
    for (size_t j = 0; j < rows; j++) {
      f += trans(i,j) * data[j];
    }
    for (size_t j = 0; j < rank; j++) {
      if (j != c) {
	PolynomialFp<R,S> x_j(this->GF, j);
	f -= mat(i,j) * x_j;
      }
    }
    eval_list[c] = f;
  }

  // The matrix whose rows parameterize all isotropic subspaces.
  for (size_t row = 0; row < p_isotropic_param->nrows(); row++)
    for (size_t col = 0; col < p_isotropic_param->ncols(); col++)
      (*p_isotropic_param)(row,col) =
	(*p_isotropic_param)(row,col).evaluate(eval_list);

#ifdef DEBUG
  // Verify that we didn't screw up somewhere along the line.
  for (size_t i = 0; i < this->k; i++)
    for (size_t j = 0; j < this->k; j++) {
      std::vector<PolynomialFp<R, S> > vec;
      for (size_t r = 0; r < n; r++)
	vec.push_back((i == j) ? (*p_isotropic_param)(i,r) :
		      (*p_isotropic_param)(i,r) + (*p_isotropic_param)(j,r));
      PolynomialFp<R, S> f = p_q_std->evaluate(vec);
      assert(f == zero);
    }
#endif

  // Determine the free variables.

  for (size_t i = 0; i < rank; i++) {
    bool appears = false;
    for (size_t row = 0; row < p_isotropic_param->nrows(); row++)
      for (size_t col = 0; col < p_isotropic_param->ncols(); col++)
	if ((*p_isotropic_param)(row,col).degree(i) > 0) {
	  appears = true;
	  break;
	}
    if (!appears)
      remove.push_back(i);
  }

  for (size_t i = 0; i < rank; i++) {
    std::vector<size_t>::const_iterator it;
    it = std::find(remove.begin(), remove.end(), i);
    if (it == remove.end()) {
      this->free_vars.push_back(i);
      this->params.push_back(zero);
    }
  }
  return;
}
