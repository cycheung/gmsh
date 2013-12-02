#include <complex>
#include <sstream>
#include "SolverMatrix.h"

using namespace std;

template<>
string SolverMatrix<double>::toMatlab(string matrixName) const{
  // Init
  stringstream stream;
  list<pair<size_t, double> >::iterator it;
  list<pair<size_t, double> >::iterator end;

  // Common part
  stream << matlabCommon(matrixName);

  // Values
  stream << "[";
  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << std::scientific << it->second << ", ";
  }
  stream << "], ";

  // Number of rows and columns
  stream << nRow << ", " << nCol << ")";

  // Return
  return stream.str();
}

template<>
string SolverMatrix<std::complex<double> >::toMatlab(string matrixName) const{
  // Init
  stringstream stream;
  list<pair<size_t, std::complex<double> > >::iterator it;
  list<pair<size_t, std::complex<double> > >::iterator end;

  // Common part
  stream << matlabCommon(matrixName);

  // Values
  stream << "[";
  for(size_t i = 0; i < nRow; i++){
    it  = data[i].begin();
    end = data[i].end();

    for(; it != end; it++)
      stream << std::scientific
             << "(" << it->second.real() << " + i * "<< it->second.imag() << ")"
             << ", ";
  }
  stream << "], ";

  // Number of rows and columns
  stream << nRow << ", " << nCol << ")";

  // Return
  return stream.str();
}
