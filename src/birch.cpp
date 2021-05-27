#include <set>

#include "birch.h"
#include "Genus.h"
#include "IsometrySequence.h"

int main(int argc, char **argv)
{
    std::cout<< "In main of birch.cpp" << std::endl;
    std::vector<Z64_PrimeSymbol> symbols_64;
    Z64_PrimeSymbol p_64;
    std::vector<Z_PrimeSymbol> symbols;
    Z_PrimeSymbol p;
    
    Z64_QuadForm<3>::SymVec coeffs_64 = {2,1,2,1,1,2};
    Z_QuadForm<3>::SymVec coeffs = {2,1,2,1,1,2};

    Z64_QuadForm<3> q0_64(coeffs_64);
    const Z64_SquareMatrix<3> & B_64 = q0_64.bilinear_form();

    std::cout << "B_64 = " << B_64 << std::endl;
    
    Z_QuadForm<3> q0(coeffs);
    const Z_SquareMatrix<3> & B = q0.bilinear_form();

    std::cout << "B = " << B << std::endl;

    std::vector<std::vector<Z64_QuadForm<5> > >
      vec_64 = Z64_QuadForm<5>::get_quinary_forms(61);
    
    std::vector<std::vector<Z_QuadForm<5> > >
      vec = Z_QuadForm<5>::get_quinary_forms(61);

    std::set<Z64> F_64;
    std::set<std::pair<Z64, int> > F_ext_64;
    
    std::set<Z> F;
    std::set<std::pair<Z, int> > F_ext;

    size_t I;

    for (std::vector<Z64_QuadForm<5> > genus : vec_64)
      {
	for (Z64_QuadForm<5> q : genus)
	  {
	    std::cout << q << std::endl;
	    std::cout << q.discriminant() << std::endl;
	    Z64 det = q.invariants(F_64,I);
	    std::cout << "det = " << det << std::endl;
	    std::cout<< std::endl;
	    for (Z64 f : F_64)
	      std::cout << f << " ";
	    std::cout<< std::endl << I << std::endl << std::endl;
	    det = q.invariants(F_ext_64,I);
	    for (std::pair<Z64,int> f : F_ext_64)
	      std::cout << "Hasse(" << f.first << ") =  " << f.second << " ";
	    std::cout<< std::endl;
	  }
      }
    
    for (std::vector<Z_QuadForm<5> > genus : vec)
      {
	for (Z_QuadForm<5> q : genus)
	  {
	    std::cout << q << std::endl;
	    std::cout << q.discriminant() << std::endl;
	    Z det = q.invariants(F,I);
	    std::cout << "det = " << det << std::endl;
	    std::cout<< std::endl;
	    for (Z f : F)
	      std::cout << f << " ";
	    std::cout<< std::endl << I << std::endl << std::endl;
	    det = q.invariants(F_ext,I);
	    for (std::pair<Z,int> f : F_ext)
	      std::cout << "Hasse(" << f.first << ") =  " << f.second << " ";
	    std::cout<< std::endl;
	  }
      }

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

    genus1->hecke_matrix_dense(8191);
    // Here we overflow because 8191^3 is very nearly 32 bits
    // How did Jeff get over that?
    //     genus2->hecke_matrix_dense(8191);

    return EXIT_SUCCESS;
}
