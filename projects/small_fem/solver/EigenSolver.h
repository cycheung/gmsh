// Gmsh - Copyright (C) 1997-2012 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to <gmsh@geuz.org>.

#ifndef _EIGENSOLVER_H_
#define _EIGENSOLVER_H_

#include <string>
#include <complex>
#include "GmshConfig.h"
#include "GmshMessage.h"

#include "linearSystemPETSc.h"
#include <slepc.h>

class EigenSolver{
 private:
  linearSystemPETSc<double> *_A, *_B;
  bool _hermitian;
  std::vector<std::complex<double> > _eigenValues;
  std::vector<std::vector<std::complex<double> > > _eigenVectors;
  void _try(int ierr) const { CHKERRABORT(PETSC_COMM_WORLD, ierr); }

 public:
  EigenSolver(linearSystemPETSc<double> *A, linearSystemPETSc<double>* B = NULL,
              bool hermitian=true);
  bool solve(int numEigenValues=0, std::string which="");
  int getNumEigenValues(){ return _eigenValues.size(); }
  std::complex<double> getEigenValue(int num){ return _eigenValues[num]; }
  std::vector<std::complex<double> > &getEigenVector(int num){ return _eigenVectors[num]; }
  void clear()
  {
    _eigenValues.clear();
    _eigenVectors.clear();
  };
  std::complex<double> getEigenVectorComp(int num, int com)
  {
    return _eigenVectors[num][com];
  };
};

#endif
