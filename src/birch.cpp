#include <set>

#include "birch.h"
#include "Genus.h"
#include "IsometrySequence.h"

int main(int argc, char **argv)
{
#ifdef DEBUG_LEVEL_FULL
    std::cerr<< "In main of birch.cpp" << std::endl;
#endif
    
    std::vector<Z64_PrimeSymbol> symbols_64;
    Z64_PrimeSymbol p_64;
    std::vector<Z_PrimeSymbol> symbols;
    Z_PrimeSymbol p;
    
    Z64_QuadForm<3>::SymVec coeffs_64 = {2,1,2,1,1,2};
    Z_QuadForm<3>::SymVec coeffs = {2,1,2,1,1,2};

    Z64_QuadForm<3> q0_64(coeffs_64);
    
#ifdef DEBUG_LEVEL_FULL
    const Z64_SquareMatrix<3> & B_64 = q0_64.bilinear_form();
    std::cerr << "B_64 = " << B_64 << std::endl;
#endif
    
    Z_QuadForm<3> q0(coeffs);
    
#ifdef DEBUG_LEVEL_FULL
    const Z_SquareMatrix<3> & B = q0.bilinear_form();
    std::cerr << "B = " << B << std::endl;
#endif
    
    std::vector<std::vector<Z64_QuadForm<5> > >
      vec_64 = Z64_QuadForm<5>::get_quinary_forms(61);
    
    std::vector<std::vector<Z_QuadForm<5> > >
      vec = Z_QuadForm<5>::get_quinary_forms(61);

    std::set<Z64> F_64;
    std::set<std::pair<Z64, int> > F_ext_64;
    
    std::set<Z> F;
    std::set<std::pair<Z, int> > F_ext;
    
#ifdef DEBUG_LEVEL_FULL
    size_t I;

    for (std::vector<Z64_QuadForm<5> > genus : vec_64)
      {
	for (Z64_QuadForm<5> q : genus)
	  {
	    std::cerr << q << std::endl;
	    std::cerr << q.discriminant() << std::endl;
	    Z64 det = q.invariants(F_64,I);
	    std::cerr << "det = " << det << std::endl;
	    std::cerr<< std::endl;
	    for (Z64 f : F_64)
	      std::cerr << f << " ";
	    std::cerr<< std::endl << I << std::endl << std::endl;
	    det = q.invariants(F_ext_64,I);
	    for (std::pair<Z64,int> f : F_ext_64)
	      std::cerr << "Hasse(" << f.first << ") =  " << f.second << " ";
	    std::cerr << std::endl;
	  }
      }
    
    for (std::vector<Z_QuadForm<5> > genus : vec)
      {

	for (Z_QuadForm<5> q : genus)
	  {
	    std::cerr << q << std::endl;
	    std::cerr << q.discriminant() << std::endl;
	    Z det = q.invariants(F,I);
	    std::cerr << "det = " << det << std::endl;
	    std::cerr<< std::endl;
	    for (Z f : F)
	      std::cerr << f << " ";
	    std::cerr<< std::endl << I << std::endl << std::endl;
	    det = q.invariants(F_ext,I);
	    for (std::pair<Z,int> f : F_ext)
	      std::cerr << "Hasse(" << f.first << ") =  " << f.second << " ";
	    std::cerr << std::endl;
	  }

      }
#endif
    
    p_64.p = 61;
    p_64.power = 1;
    p_64.ramified = true;
    symbols_64.push_back(p_64);
    
    p.p = 61;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    Z64_QuadForm<5> q5_64 = vec_64[0][0];
    Z64_Genus<5> genus5_64(q5_64, symbols_64);
    
    Z_QuadForm<5> q5 = vec[0][0];
    Z_Genus<5> genus5(q5, symbols);
    
    // we start by just considering the single prime 61
    /* 
    p.p = 11;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    p.p = 13;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    p.p = 17;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    p.p = 19;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);
  
    p.p = 23;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);
    */
    Z_Isometry<3> s;
    Z_QuadForm<3> q = Z_QuadForm<3>::get_quad_form(symbols);

    Z_Genus<3> genus1(q, symbols);
    std::shared_ptr<Z64_Genus<3> > genus2 = std::make_shared<Z64_Genus<3> >(genus1);

    //  std::map<Z64, std::vector<std::vector<int> > > T2 =
    //  genus2->hecke_matrix_sparse(2);
    // std::map<Z64, std::vector<std::vector<int> > >::const_iterator i;
    std::vector<Z64> primes = {2,3,5,7,11,13,17,19};
    for(size_t j = 0; j < primes.size(); j++) {
      std::map<Z64, std::vector<int> > T =
	genus2->hecke_matrix_dense(primes[j]);
      std::map<Z64, std::vector<int> >::const_iterator i;
      for (i = T.begin(); i != T.end(); i++) {
	std::cout << " with spinor " << i->first << std::endl;
	std::cout << " T_" << primes[j] << " = " << i->second << std::endl;
      }
    }

    
    for(size_t j = 0; j < primes.size(); j++) {
      // Here we have an overflow already at 11
      // std::map<Z64, std::vector<int> > T =
      //	genus5_64.hecke_matrix_dense(primes[j]);
      //      std::map<Z64, std::vector<int> >::const_iterator i;
      std::map<Z, std::vector<int> > T =
	genus5.hecke_matrix_dense(primes[j]);
      std::map<Z, std::vector<int> >::const_iterator i;
      for (i = T.begin(); i != T.end(); i++) {
	std::cout << " with spinor " << i->first << std::endl;
	std::cout << " T_" << primes[j] << " = " << i->second << std::endl;
      }
    }
    // genus1.hecke_matrix_dense(8191);
    
    // Here we overflow because 8191^3 is very nearly 32 bits
    // How did Jeff get over that?
    //     genus2->hecke_matrix_dense(8191);

    return EXIT_SUCCESS;
}
