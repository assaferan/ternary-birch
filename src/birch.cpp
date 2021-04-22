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
    const Z_QuadForm::RMat & B = q0.getBilinearForm();

    for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
          std::cout << B[i][j] << " ";
	std::cout << std::endl;
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
