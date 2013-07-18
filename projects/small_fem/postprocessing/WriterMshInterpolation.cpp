#include "WriterMsh.h"

using namespace std;

void WriterMsh::writeInterpolationScheme(void) const{
  // Some Temp Value
  const fullMatrix<double>& coef = lBasis->getCoefficient();
  const fullMatrix<double>& mono = lBasis->getMonomial();

  const size_t nRowCoef = coef.size1();
  const size_t nColCoef = coef.size2();

  const size_t nRowMono = mono.size1();
  const size_t nColMono = mono.size2();

  // Up to now, we suppose *ONE* topology
  *out << "$InterpolationScheme"     << endl
       << "\"interpolation scheme\"" << endl
       << "1"                        << endl
       << (*element)[0]->getType()   << endl

    // 2 Matrices: Coefficients and Monomials
       << "2"                        << endl;

  // Coefficients Matrix
  *out << nRowCoef << " "
       << nColCoef << endl;

  for(size_t i = 0; i < nRowCoef; i++){
    for(size_t j = 0; j < nColCoef; j++){
      *out << coef(i, j);

      if(j < nColCoef - 1)
        *out << " ";

      else
        *out << endl;
    }
  }

  // Monomials Matrix
  *out << nRowMono << " "
       << nColMono << endl;

  for(size_t i = 0; i < nRowMono; i++){
    for(size_t j = 0; j < nColMono; j++){
      *out << mono(i, j);

      if(j < nColMono - 1)
        *out << " ";

      else
        *out << endl;
    }
  }

  // End
  *out << "$EndInterpolationScheme" << endl;
}
