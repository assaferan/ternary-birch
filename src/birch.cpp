#include <set>

#include "birch.h"
#include "Genus.h"
#include "IsometrySequence.h"

int main(int argc, char **argv)
{
    std::cout<< "In main of birch.cpp" << std::endl;
    std::vector<Z_PrimeSymbol> symbols;
    Z_PrimeSymbol p;

    Z_QuadForm<3>::RVec coeffs = {2,1,2,1,1,2};

    Z_QuadForm<3> q0(coeffs);
    const Z_QuadForm<3>::RMat & B = q0.bilinear_form();

    for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
          std::cout << B[i][j] << " ";
	std::cout << std::endl;
      }

    std::vector<std::vector<Z_QuadForm<5> > >
      vec = Z_QuadForm<5>::get_quinary_forms(61);

    std::set<Z> F;
    std::set<std::pair<Z, int> > F_ext;

    /*
    size_t I;
    
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
    */
    p.p = 61;
    p.power = 1;
    p.ramified = true;
    symbols.push_back(p);

    Z_QuadForm<5> q5 = vec[0][0];
    Z_Genus<5> genus5(q5, symbols);
    
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

    Z_Isometry<3> s;
    Z_QuadForm<3> q = Z_QuadForm<3>::get_quad_form(symbols);

    Z_Genus<3> genus1(q, symbols);
    std::shared_ptr<Z64_Genus<3> > genus2 = std::make_shared<Z64_Genus<3> >(genus1);

    genus2->hecke_matrix_dense(8191);

    return EXIT_SUCCESS;
}
