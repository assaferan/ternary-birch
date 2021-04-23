#include "birch.h"
#include "Genus.h"
#include "IsometrySequence.h"

int main(int argc, char **argv)
{
    std::cout<< "In main of birch.cpp" << std::endl;
    std::vector<Z_PrimeSymbol> symbols;
    Z_PrimeSymbol p;

    Z_QuadForm::RVec coeffs = {2,1,2,1,1,2};

    Z_QuadForm q0(coeffs);
    const Z_QuadForm::RMat & B = q0.bilinear_form();

    for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
          std::cout << B[i][j] << " ";
	std::cout << std::endl;
      }

    std::vector<std::vector<QuadForm<Z,5> > >
      vec = QuadForm<Z,5>::get_quinary_forms(61);

    const typename QuadForm<Z,5>::RDiag & D;
    std::set<Z> F;
    size_t I;
    
    for (std::vector<QuadForm<Z,5> > genus : vec)
      {
	for (QuadForm<Z,5> q : genus)
	  {
	    std::cout << q << std::endl;
	    std::cout << q.discriminant() << std::endl;
	    q.invariants(D,F,I);
	    //	    const typename QuadForm<Z,5>::RDiag & d = q.orthogonalize_gram();
	    for (size_t j = 0; j < 5; j++)
	      std::cout << D[j] << " ";
	    std::cout<< std::endl;
	    for (R f : F)
	      std::cout << f << " ";
	    std::cout<< std::endl << I << std::endl << std::endl;
	    
	  }
      }
    
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

    Z_Isometry s;
    Z_QuadForm q = Z_QuadForm::get_quad_form(symbols);

    Z_Genus genus1(q, symbols);
    std::shared_ptr<Z64_Genus> genus2 = std::make_shared<Z64_Genus>(genus1);

    genus2->hecke_matrix_dense(8191);

    return EXIT_SUCCESS;
}
