#include "WriterMsh.h"

using namespace std;

void WriterMsh::writeInterpolationScheme(void) const{
  // Some Temp Value
  const fullMatrix<double>& coef = lBasis->getCoefficient();
  const fullMatrix<double>& mono = lBasis->getMonomial();

  const unsigned int nRowCoef = coef.size1();
  const unsigned int nColCoef = coef.size2();

  const unsigned int nRowMono = mono.size1();
  const unsigned int nColMono = mono.size2();

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

  for(unsigned int i = 0; i < nRowCoef; i++){
    for(unsigned int j = 0; j < nColCoef; j++){
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

  for(unsigned int i = 0; i < nRowMono; i++){
    for(unsigned int j = 0; j < nColMono; j++){
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
