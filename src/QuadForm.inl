// Implementation of templated functions from QuadForm.h
#include "Math.h"

// c-tors
template<typename R, size_t n>
QuadForm<R, n>::QuadForm(const RVec& coeffs)
{
  size_t idx = 0;
  for (size_t row = 0; row < n; row++)
    {
      for (size_t col = 0; col < row; col++)
	{	
	  this->B_[row][col] = coeffs[idx];
	  this->B_[col][row] = coeffs[idx++];
	}
      this->B_[row][row] = coeffs[idx++];
    }
}

template<typename R, size_t n>
R QuadForm<R, n>::discriminant(void) const
{
  if (n == 3) 
        return this->a_ * (4 * this->b_ * this->c_ - this->f_ * this->f_) -
            this->b_ * this->g_ * this->g_ +
            this->h_ * (this->f_ * this->g_ - this->c_ * this->h_);
  else
  {
  // Instead of the previous ad-hoc method, we use Bareiss algorithm
  // to compute the determinant.
  // TODO - can do Cholesky, will be faster
  // !! TODO - Leibniz should be better when n <= 5
    R M[n+1][n+1];
    M[0][0] = 1;
    // init
    for (size_t row = 0; row < n; row++)
      for (size_t col = 0; col < n; col++)
        M[row+1][col+1] = this->B_[row][col];
    for (size_t k = 1; k < n; k++)
      for (size_t i = k+1; i <= n; i++)
        for (size_t j = k+1; j <= n; j++)
          M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j])/M[k-1][k-1];
    if (n % 2 == 0)	  
      return M[n][n];
    else
      return M[n][n]/2;
  }
}

template<typename R, size_t n>
std::vector<R> QuadForm<R, n>::orthogonalize_gram() const
{
  std::vector<R> D(n);
  typename QuadForm<R,n>::RMat L;
  R prod_diag = 1;
  R d, inner_sum;
  // This works but inefficiently - for some reason we get O(n^4) operations.
  // !! TODO - check it out later
  // Oh I see - we should do the L update in two passes...
  for (size_t i = 0; i < n; i++)
    {
      L[i][i] = prod_diag;
      d = prod_diag;
      for (size_t j = 0; j < i; j++)
	{
	  L[i][j] = 0;
	  for (size_t k = 0; k < i; k++)
	    {
	      inner_sum = 0;
	      for (size_t r = 0; r <= k; r++)
		inner_sum += L[k][r]*(this->B_[i][r])*L[k][j];
	      inner_sum *= -L[i][i] / D[k];
	      L[i][j] += inner_sum;
	    }
	  d = gcd(d, L[i][j]);
	}
      for (size_t j = 0; j <= i; j++)
	L[i][j] /= d;
      D[i] = 0;
      for (size_t j = 0; j <= i; j++)
	for (size_t k = 0; k <= i; k++)
	  D[i] += L[i][j]*(this->B_[j][k])*L[i][k];
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
int QuadForm<R,n>::Hasse(const std::vector<R> & D, const R & p)
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
R QuadForm<R, n>::invariants(std::set<R> & F, size_t& I) const
{
  std::vector<R> D = this->orthogonalize_gram();
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
R QuadForm<R, n>::invariants(std::set<std::pair<R, int> > & F, size_t& I) const
{
  std::vector<R> D = this->orthogonalize_gram();
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
    F.insert(std::make_pair(p, Hasse(D,p)));

  R prod = 1;
  for (size_t i = 0; i < n; i++)
    prod *= D[i];
  
  return prod;
}

template<typename R, size_t n>
R QuadForm<R,n>::inner_product(const typename QuadForm<R,n>::RMat & F,
			       const typename QuadForm<R,n>::RMat & S,
			       size_t idx1, size_t idx2)
{
  R ans = 0;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      ans += S[idx1][i] * F[i][j] * S[idx2][j];
  return ans;
}

template<typename R, size_t n>
typename QuadForm<R, n>::jordan_data
QuadForm<R, n>::jordan_decomposition(const R & p) const
{
  bool even = (p == 2);
  QuadForm<R, n>::RMat S, G;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      S[i][j] = (i == j) ? 1 : 0;
  size_t k = 0;
  // virtually infinity
  size_t old_val = 0xffffffff;
  std::vector<size_t> blocks;
  jordan_data jordan;
  while (k < n)
    {
      // G = SFS^t
     for (size_t i = 0; i < n; i++)
       for (size_t j = 0; j < n; j++)
	   G[i][j] = inner_product(this->B_, S, i, j);

     size_t ii = k;
     size_t m = Math<R>::valuation(G[k][k], p);
     
     for (size_t i = k+1; i < n; i++)
       {
	 size_t val = Math<R>::valuation(G[i][i], p);
	 if (val < m)
	   {
	     m = val;
	     ii = i;
	   }
       }
     std::pair<size_t, size_t> i_pair = std::make_pair(ii, ii);
     for (size_t i = k; i < n; i++)
       for (size_t j = i+1; j < n; j++)
	 {
	   size_t tmp = Math<R>::valuation(G[i][j], p);
	   if (tmp < m)
	     {
	       m = tmp;
	       i_pair.first = i;
	       i_pair.second = j;
	     }
	 }
     if (m != old_val)
       {
	 blocks.push_back(k);
	 old_val = m;
	 jordan.exponents.push_back(m);
       }
     if ((even) && (i_pair.first != i_pair.second))
       {
	 // swap rows
	 for (size_t i = 0; i < n; i++)
	   {
	     R tmp = S[i_pair.first][i];
	     S[i_pair.first][i] = S[k][i];
	     S[k][i] = tmp;
	   }
	 // swap rows
	 for (size_t i = 0; i < n; i++)
	   {
	     R tmp = S[i_pair.second][i];
	     S[i_pair.second][i] = S[k+1][i];
	     S[k+1][i] = tmp;
	   }
	 // T12 = S[k]*F*S[k+1]^t
	 R T12 = inner_product(this->B_, S, k, k+1);

	 // multiply S[k] by p^val(T12,p)/T12
	 // Check whether we have to change to rational here
	 for (size_t i = 0; i < n; i++)
	   S[k][i] *= (1 << Math<R>::valuation(T12, p)) / T12;
	 R T11 = inner_product(this->B_, S, k, k);
	 R T22 = inner_product(this->B_, S, k+1, k+1);
	 T12 = inner_product(this->B_, S, k, k+1);
	 R d = T11*T22-T12*T12;
	 for (size_t l = k+2; l < n; l++)
	   {
	     R tl = T12*inner_product(this->B_,S,k+1,l) -
	       T22*inner_product(this->B_,S,k,l);
	     R ul = T12*inner_product(this->B_,S,k,l) -
	       T11*inner_product(this->B_,S,k+1,l);
	     for (size_t i = 0; i < n; i++)
	       S[l][i] += (tl/d)*S[k][i] + (ul/d)*S[k+1][i];
	   }
	 k += 2;
       }
     else
       {
	 if (i_pair.first == i_pair.second)
	   // swap rows
	   for (size_t i = 0; i < n; i++)
	     {
	       R tmp = S[i_pair.first][i];
	       S[i_pair.first][i] = S[k][i];
	       S[k][i] = tmp;
	     }
	 else
	   {
	     for (size_t i = 0; i < n; i++)
	       S[i_pair.first][i] += S[i_pair.second][i];
	     // swap rows
	     for (size_t i = 0; i < n; i++)
	     {
	       R tmp = S[i_pair.first][i];
	       S[i_pair.first][i] = S[k][i];
	       S[k][i] = tmp;
	     }
	   }
	 R nrm = inner_product(this->B_, S, k, k);
	 R X[n];
	 for (size_t i = 0; i < n; i++)
	   X[i] = inner_product(this->B_, S, k, i);
	 for (size_t l = k+1; l < n; l++)
	     for (size_t i = 0; i < n; i++)
	       S[l][i] -= X[l]/nrm * S[k][i];
	 k += 1;
       }
    }
  blocks.push_back(n+1);

  for (size_t i = 0; i < blocks.size()-1; i++) {
    size_t nrows = blocks[i+1]-blocks[i];
    std::vector<R> data(nrows*n);
    size_t idx = 0;
    for (size_t row = 0; row < nrows; row++)
      for (size_t col = 0; col < n; col++)
	data[idx++] = S[blocks[i]+row][col];
    Matrix<R> mat(data, nrows, n);
    jordan.matrices.push_back(mat);
  }
  for (Matrix<R> m  : jordan.matrices) {
    Matrix<R> F(this->B_);
    jordan.grams.push_back(m*F*m.transpose());
  }
  return jordan;
}

template<typename R, size_t n>
int QuadForm<R, n>::border(const QuadForm<R, n>& q, int m)
{	     
  switch (m)
    {
    case 1:
      return (q.a() == q.h()) && (q.g() == q.f()*2);
    case 2:
      return (q.a() == q.g()) && (q.h() == q.f()*2);
    case 3:
      return (q.b() == q.f()) && (q.h() == q.g()*2);
    case 4:
      return (q.a() == -q.h());
    case 5:
      return (q.a() == -q.g());
    case 6:
      return (q.b() == -q.f());
    case 7:
      return (q.a() + q.b() + q.f() + q.g() + q.h() == 0) &&
	(q.a()*2 + q.g()*2 + q.h() == 0);
    case 8:
      return (q.a() == q.b()) && (q.f() == q.g());
    case 9:
      return (q.b() == q.c()) && (q.g() == q.h());
    case 10:
      return (q.f() == q.g()) && (q.f() == 0);
    case 11:
      return (q.f() == q.h()) && (q.f() == 0);
    case 12:
      return (q.g() == q.h()) && (q.g() == 0);
    case 13:
      return (q.f() == q.g()) && (q.g() == q.h()) &&
	(q.h() == q.a());
    case 14:
      return (q.a() == q.g()) && (q.a() == q.h());
    case 15:
      return (q.a() == q.b()) &&
	(q.a() + q.b() + q.f() + q.g() + q.h() == 0);
    case 16:
      return (q.a() == q.b()) && (q.b() == q.c()) &&
	(q.a() + q.b() + q.f() + q.g() + q.h() == 0);
    default:
      return 0;
    }
}

template<typename R, size_t n>
int QuadForm<R, n>::num_automorphisms(const QuadForm<R, n>& q)
{
  if (border(q, 1))
    {
      if (border(q, 2))
	{
	  if (border(q, 14))
	    {
	      if (border(q, 9))
		return 16;
	      else
		return 8;
	    }
	}
      else
	return 4;
    }

  if (border(q, 2))
    return 4;

  if (border(q, 3))
    return 4;

  if (border(q, 4))
    {
      if (border(q, 10))
	{
	  if (border(q, 8))
	    return 24;
	  else
	    return 8;
	}
      else
	return 4;
    }

  if (border(q, 5))
    {
      if (border(q, 6))
	{
	  if (border(q, 7))
	    {
	      if (border(q, 8))
		{
		  if (border(q, 15))
		    return 16;
		}
	      else
		return 8;
	    }
	}
      else if (border(q, 11))
	return 8;
      else
	return 4;
    }

  if (border(q, 6))
    {
      if (border(q, 12))
	{
	  if (border(q, 9))
	    return 24;
	  else
	    return 8;
	}
      else
	return 4;
    }

  if (border(q, 7))
    {
      if (border(q, 8) && border(q, 15))
	{
	  if (border(q, 16))
	    {
	      if (border(q, 9))
		return 48;
	      else
		return 16;
	    }
	  else
	    return 8;
	}
      else if (border(q, 9))
	return 12;
      else
	return 4;
    }

  if (border(q, 8))
    {
      if (border(q, 9))
	{
	  if (border(q, 10) && border(q, 11) && border(q, 12))
	    return 48;
	  else if (border(q, 13) && border(q, 14))
	    return 48;
	  else
	    return 12;
	}
      else if (border(q, 10))
	{
	  if (border(q, 11) && border(q, 12))
	    return 16;
	  else
	    return 8;
	}
      else if (border(q, 14))
	return 12;
      else
	return 4;
    }

  if (border(q, 9))
    {
      if (border(q, 12))
	{
	  if (border(q, 10) && border(q, 11))
	    return 16;
	  else
	    return 8;
	}
      else if (border(q, 14))
	{
	  if (border(q, 13))
	    return 8;
	  else
	    return 8;
	}
      else if (border(q, 15))
	return 16;
      else
	return 4;
    }

  if (border(q, 10))
    {
      if (border(q, 11) && border(q, 12))
	return 8;
      else
	return 4;
    }

  if (border(q, 11))
    return 4;

  if (border(q, 12))
    return 4;

  if (border(q, 13) && border(q, 14))
    return 4;

  if (border(q, 14))
    return 4;

  if (border(q, 15))
    {
      if (border(q, 16))
	return 8;
      else
	return 4;
    }

  return 2;
}

template<typename R, size_t n>
const std::vector<Isometry<R,n>>&
QuadForm<R,n>::proper_automorphisms(const QuadForm<R, n>& q)
{
  if (border(q, 1))
    {
      if (border(q, 2))
	{
	  if (border(q, 14))
	    {
	      if (border(q, 9))
		{
		  return Isometry<R,n>::automorphisms[0];
		}
	      else
		{
		  return Isometry<R,n>::automorphisms[1];
		}
	    }
	}
      else
	{
	  return Isometry<R,n>::automorphisms[2];
	}
    }

  if (border(q, 2))
    {
      return Isometry<R,n>::automorphisms[3];
    }

  if (border(q, 3))
    {
      return Isometry<R,n>::automorphisms[4];
    }

  if (border(q, 4))
    {
      if (border(q, 10))
	{
	  if (border(q, 8))
	    {
	      return Isometry<R,n>::automorphisms[5];
	    }
	  else
	    {
	      return Isometry<R,n>::automorphisms[6];
	    }
	}
      else
	{
	  return Isometry<R,n>::automorphisms[7];
	}
    }

  if (border(q, 5))
    {
      if (border(q, 6))
	{
	  if (border(q, 7))
	    {
	      if (border(q, 8))
		{
		  if (border(q, 15))
		    {
		      return Isometry<R,n>::automorphisms[8];
		    }
		}
	      else
		{
		  return Isometry<R,n>::automorphisms[9];
		}
	    }
	}
      else if (border(q, 11))
	{
	  return Isometry<R,n>::automorphisms[10];
	}
      else
	{
	  return Isometry<R,n>::automorphisms[11];
	}
    }

  if (border(q, 6))
    {
      if (border(q, 12))
	{
	  if (border(q, 9))
	    {
	      return Isometry<R,n>::automorphisms[12];
	    }
	  else
	    {
	      return Isometry<R,n>::automorphisms[13];
	    }
	}
      else
	{
	  return Isometry<R,n>::automorphisms[14];
	}
    }

  if (border(q, 7))
    {
      if (border(q, 8) && border(q, 15))
	{
	  if (border(q, 16))
	    {
	      if (border(q, 9))
		{
		  return Isometry<R,n>::automorphisms[15];
		}
	      else
		{
		  return Isometry<R,n>::automorphisms[16];
		}
	    }
	  else
	    {
	      return Isometry<R,n>::automorphisms[17];
	    }
	}
      else if (border(q, 9))
	{
	  return Isometry<R,n>::automorphisms[18];
	}
      else
	{
	  return Isometry<R,n>::automorphisms[19];
	}
    }

  if (border(q, 8))
    {
      if (border(q, 9))
	{
	  if (border(q, 10) && border(q, 11) && border(q, 12))
	    {
	      return Isometry<R,n>::automorphisms[20];
	    }
	  else if (border(q, 13) && border(q, 14))
	    {
	      return Isometry<R,n>::automorphisms[21];
	    }
	  else
	    {
	      return Isometry<R,n>::automorphisms[22];
	    }
	}
      else if (border(q, 10))
	{
	  if (border(q, 11) && border(q, 12))
	    {
	      return Isometry<R,n>::automorphisms[23];
	    }
	  else
	    {
	      return Isometry<R,n>::automorphisms[24];
	    }
	}
      else if (border(q, 14))
	{
	  return Isometry<R,n>::automorphisms[25];
	}
      else
	{
	  return Isometry<R,n>::automorphisms[26];
	}
    }

  if (border(q, 9))
    {
      if (border(q, 12))
	{
	  if (border(q, 10) && border(q, 11))
	    {
	      return Isometry<R,n>::automorphisms[27];
	    }
	  else
	    {
	      return Isometry<R,n>::automorphisms[28];
	    }
	}
      else if (border(q, 14))
	{
	  if (border(q, 13))
	    {
	      return Isometry<R,n>::automorphisms[29];
	    }
	  else
	    {
	      return Isometry<R,n>::automorphisms[30];
	    }
	}
      else if (border(q, 15))
	{
	  return Isometry<R,n>::automorphisms[31];
	}
      else
	{
	  return Isometry<R,n>::automorphisms[32];
	}
    }

  if (border(q, 10))
    {
      if (border(q, 11) && border(q, 12))
	{
	  return Isometry<R,n>::automorphisms[33];
	}
      else
	{
	  return Isometry<R,n>::automorphisms[34];
	}
    }

  if (border(q, 11))
    {
      return Isometry<R,n>::automorphisms[35];
    }

  if (border(q, 12))
    {
      return Isometry<R,n>::automorphisms[36];
    }

  if (border(q, 13) && border(q, 14))
    {
      return Isometry<R,n>::automorphisms[37];
    }

  if (border(q, 14))
    {
      return Isometry<R,n>::automorphisms[38];
    }

  if (border(q, 15))
    {
      if (border(q, 16))
	{
	  return Isometry<R,n>::automorphisms[39];
	}
      else
	{
	  return Isometry<R,n>::automorphisms[40];
	}
    }

  return Isometry<R,n>::automorphisms[41];
}

template<typename R, size_t n>
QuadForm<R,n> QuadForm<R,n>::reduce(const QuadForm<R,n>& q, Isometry<R,n>& s)
{
  R a = q.a_;
  R b = q.b_;
  R c = q.c_;
  R f = q.f_;
  R g = q.g_;
  R h = q.h_;

  int flag = 1;
  while (flag)
    {
      R t = a + b + f + g + h;
      if (t < 0)
	{
	  s.A101011001();
	  c += t;
	  f += (h + b + b);
	  g += (h + a + a);
	}

      if (a >= h)
	{
	  t = birch_util::dumb_div<R>(a-h, a+a);
	}
      else
	{
	  t = -birch_util::dumb_div<R>(a+h-1, a+a);
	}
      if (t != 0)
	{
	  s.A1t0010001(t);
	  R temp = a * t;
	  h += temp;
	  b += h * t;
	  f += g * t;
	  h += temp;
	}

      if (b >= f)
	{
	  t = birch_util::dumb_div<R>(b-f, b+b);
	}
      else
	{
	  t = -birch_util::dumb_div<R>(b+f-1, b+b);
	}
      if (t != 0)
	{
	  s.A10001t001(t);
	  R temp = b * t;
	  f += temp;
	  c += f * t;
	  g += h * t;
	  f += temp;
	}

      if (a >= g)
	{
	  t = birch_util::dumb_div<R>(a-g, a+a);
	}
      else
	{
	  t = -birch_util::dumb_div<R>(a+g-1, a+a);
	}
      if (t != 0)
	{
	  s.A10t010001(t);
	  R temp = a * t;
	  g += temp;
	  c += g * t;
	  f += h * t;
	  g += temp;
	}

      if (a > b || (a == b && abs(f) > abs(g)))
	{
	  s.A0n0n0000n();
	  t = a; a = b; b = t;
	  t = f; f = g; g = t;
	}

      if (b > c || (b == c && abs(g) > abs(h)))
	{
	  s.An0000n0n0();
	  t = b; b = c; c = t;
	  t = g; g = h; h = t;
	}

      if (a > b || (a == b && abs(f) > abs(g)))
	{
	  s.A0n0n0000n();
	  t = a; a = b; b = t;
	  t = f; f = g; g = t;
	}

      int fgh = (f != 0 && g != 0 && h != 0);
      if (fgh)
	{
	  if (f < 0) fgh = !fgh;
	  if (g < 0) fgh = !fgh;
	  if (h < 0) fgh = !fgh;
	}

      if (fgh)
	{
	  if (f < 0)
	    {
	      s.An00010001();
	      f = -f;
	    }

	  if (g < 0)
	    {
	      s.A1000n0001();
	      g = -g;
	    }

	  if (h < 0)
	    {
	      s.A10001000n();
	      h = -h;
	    }
	}
      else
	{
	  int s1 = f > 0;
	  int s2 = g > 0;
	  int s3 = h > 0;

	  if ((s1+s2+s3) % 2 == 1)
	    {
	      if (f == 0) s1 = 1;
	      else
		{
		  if (g == 0) s2 = 1;
		  else if (h == 0) s3 = 1;
		}
	    }

	  if (s1 == 1)
	    {
	      s.An00010001();
	      f = -f;
	    }

	  if (s2 == 1)
	    {
	      s.A1000n0001();
	      g = -g;
	    }

	  if (s3 == 1)
	    {
	      s.A10001000n();
	      h = -h;
	    }
	}

      flag = !(abs(f) <= b && abs(g) <= a &&
	       abs(h) <= a && a+b+f+g+h >= 0);
    }

  if (a + b + f + g + h == 0 &&
      a + a + g + g + h > 0)
    {
      s.An010n1001();
      c += a + b + f + g + h;
      f += h + b + b; f = -f;
      g += h + a + a; g = -g;
    }

  if (a == -h && g != 0)
    {
      s.Ann00n0001();
      f += g; f = -f;
      g = -g;
      h = -h;
    }

  if (a == -g && h != 0)
    {
      s.An0n01000n();
      f += h; f = -f;
      h = -h;
      g += (2*a);
    }

  if (b == -f && h != 0)
    {
      s.A1000nn00n();
      g += h; g = -g;
      h = -h;
      f += (2*b);
    }

  if (a == h && g > f + f)
    {
      s.Ann001000n();
      f = g - f;
    }

  if (a == g && h > f + f)
    {
      s.An0n0n0001();
      f = h - f;
    }

  if (b == f && h > g + g)
    {
      s.An000nn001();
      g = h - g;
    }

  if (a == b && abs(f) > abs(g))
    {
      s.A0n0n0000n();
      R t;
      t = a; a = b; b = t;
      t = g; g = f; f = t;
    }

  if (b == c && abs(g) > abs(h))
    {
      s.An0000n0n0();
      R t;
      t = g; g = h; h = t;
    }

  if (a == b && abs(f) > abs(g))
    {
      s.A0n0n0000n();
      R t;
      t = g; g = f; f = t;
    }

  return QuadForm<R, n>(a, b, c, f, g, h);
}

template<typename R, size_t n>
std::ostream& operator<<(std::ostream& os, const QuadForm<R,n>& q)
{
  /*
    os << "QuadForm(" << q.a_ << "," << q.b_ << "," << q.c_ << ","
    << q.f_ << "," << q.g_ << "," << q.h_ << ")";
  */
  for (size_t i = 0; i < n; i++)
    {
      for (size_t j = 0; j < n; j++)
	std::cout << q.B_[i][j] << " ";
      std::cout << std::endl;
    }
  return os;
}

template<typename R, typename S>
R QuadFormFp_3<R,S>::discriminant(void) const
{
  R res = GF->mul(this->b_, this->c_);    // bc
  res = GF->add(res, res);                // 2bc
  res = GF->add(res, res);                // 4bc
  R temp = GF->mul(this->f_, this->f_);   // ff
  res = GF->sub(res, temp);               // 4bc-ff
  res = GF->mul(this->a_, res);           // a(4bc-ff)

  temp = GF->mul(this->g_, this->g_);     // gg
  temp = GF->mul(this->b_, temp);         // bgg
  res = GF->sub(res, temp);               // a(4bc-ff)-bgg

  temp = GF->mul(this->f_, this->g_);     // fg
  R temp2 = GF->mul(this->c_, this->h_);  // ch
  temp = GF->sub(temp, temp2);            // fg-ch
  temp = GF->mul(this->h_, temp);         // h(fg-ch)
  res = GF->add(res, temp);               // a(4bc-ff)-bgg+h(fg-ch)

  return res;
}

template<typename R, typename S>
R QuadFormFp_3<R, S>::evaluate(const R& x, const R& y, const R& z) const
{
  R res = GF->mul(this->a_, x);   // ax
  R temp = GF->mul(this->g_, z);  // gz
  res = GF->add(res, temp);       // ax+gz
  temp = GF->mul(this->h_, y);    // hy
  res = GF->add(res, temp);       // ax+gz+hy
  res = GF->mul(x, res);          // x(ax+gz+hy)

  temp = GF->mul(this->b_, y);    // by
  R temp2 = GF->mul(this->f_, z); // fz
  temp = GF->add(temp, temp2);    // by+fz
  temp = GF->mul(y, temp);        // y(by+fz)
  res = GF->add(res, temp);       // x(ax+gz+hy)+y(by+fz)

  temp = GF->mul(this->c_, z);    // cz
  temp = GF->mul(temp, z);        // czz
  res = GF->add(res, temp);       // x(ax+gz+hy)+y(by+fz)+czz

  return res;
}

template<typename R, typename S, size_t n>
Vector<R,n> QuadFormFp<R, S, n>::isotropic_vector(void) const
{
  Vector<R,n> vec = {0};

  // stub - !! TODO !! - complete
  
  return vec;
}

template<typename R, typename S>
Vector3<R> QuadFormFp_3<R, S>::isotropic_vector(void) const
{
  Vector3<R> vec = {0,0,0};

  if (GF->prime() == 2) return this->isotropic_vector_p2();

  while (1)
    {
      R r = 0;
      R alpha = 0;
      R beta = 0;

      while (alpha == 0)
	{
	  r = GF->random();                   // r
	  beta = GF->mul(this->b_, r);        // br
	  alpha = GF->add(beta, this->h_);    // br+h
	  alpha = GF->mul(alpha, r);          // (br+h)r
	  alpha = GF->add(alpha, this->a_);   // (br+h)r+a = Q(1,r,0)
	  alpha = GF->mod(alpha);
	}

      R s = GF->random();

      beta = GF->add(beta, beta);         // 2br
      beta = GF->add(beta, this->h_);     // 2br+h
      beta = GF->mul(beta, s);            // (2br+h)s
      R temp = GF->mul(this->f_, r);      // fs
      beta = GF->add(beta, temp);         // (2br+h)s+fr = (dQ/dy)(s,rs,r)
      beta = GF->add(beta, this->g_);     // (2br+h)s+fr+g

      R gamma = GF->mul(this->b_, s);     // bs
      gamma = GF->add(gamma, this->f_);   // bs+f
      gamma = GF->mul(gamma, s);          // (bs+f)s
      gamma = GF->add(gamma, this->c_);   // (bs+f)s+c = Q(0,s,1)

      R disc = GF->mul(beta, beta);
      gamma = GF->mul(gamma, alpha);
      gamma = GF->add(gamma, gamma);
      gamma = GF->add(gamma, gamma);
      disc = GF->sub(disc, gamma);

      if (GF->legendre(disc) >= 0)
	{
	  R root = GF->sqrt(disc);

	  root = GF->sub(root, beta);
	  alpha = GF->add(alpha, alpha);
	  alpha = GF->mod(alpha);

	  vec.x = GF->mul(root, GF->inverse(alpha));
	  vec.y = GF->mul(r, vec.x);
	  vec.y = GF->add(vec.y, s);
	  vec.z = 1;

	  return vec;
	}
    }

  return vec;
}

// To avoid unnecessary computation, we encode each of the three 2-isotropic
  // vectors as a coordinate of the return vector. Special care must be taken
  // to obtain the actual isotropic vectors when needed.
template<typename R, typename S>
Vector3<R> QuadFormFp_3<R, S>::isotropic_vector_p2(void) const
{
  Vector3<R> vec = {0,0,0};
  R temp[3] = {0,0,0};

  int index = 0;

  if (this->c_ == 0) temp[index++] = 1;
  if (this->b_ == 0) temp[index++] = 2;
  if ((this->b_ ^ this->f_ ^ this->c_) == 0) temp[index++] = 3;
  if (this->a_ == 0) temp[index++] = 4;
  if ((this->a_ ^ this->g_ ^ this->c_) == 0) temp[index++] = 5;
  if ((this->a_ ^ this->h_ ^ this->b_) == 0) temp[index++] = 6;
  if ((this->a_ ^ this->b_ ^ this->c_ ^ this->f_ ^ this->g_ ^ this->h_) == 0) temp[index++] = 7;

  vec.x = temp[0];
  vec.y = temp[1];
  vec.z = temp[2];

  return vec;
}
