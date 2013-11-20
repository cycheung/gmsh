#ifndef _SYSTEMEIGEN_H_
#define _SYSTEMEIGEN_H_

#include <complex>
#include "SystemAbstract.h"
#include "petscmat.h"

/**
   @class SystemEigen
   @brief This class assembles an eigenvalue system

   This class assembles an eigenvalue system described by a Formulation.

   The eigenvalue problem can be generalized or not:
   @li An eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$
   @li A Generalized Eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$

   The Solver used is <a href="http://www.grycap.upv.es/slepc/">SLEPc</a>.
 */

class SystemEigen: public SystemAbstract{
 private:
  bool general;

  Mat* A;
  Mat* B;

  PetscInt nEigenValues;
  std::vector<std::complex<double> >* eigenValue;
  std::vector<fullVector<std::complex<double> > >* eigenVector;

 public:
  SystemEigen(const Formulation& formulation);
  virtual ~SystemEigen(void);

  bool isGeneral(void) const;

  size_t getEigenValuesNumber(void) const;

  const std::vector<std::complex<double> >&
    getEigenValues(void)  const;

  const std::vector<fullVector<std::complex<double> > >&
    getEigenVectors(void) const;

  void setNumberOfEigenValues(size_t nEigenValues);

  virtual void assemble(void);
  virtual void solve(void);
};


/**
   @fn SystemEigen::SystemEigen
   @param formulation A Formulation that
   gives the way to assemble the eigenvalue system

   Instantiates a new SystemEigen
   ***

   @fn SystemEigen::~SystemEigen
   Deletes this SystemEigen
   **

   @fn SystemEigen::isGeneral
   @return Returns:
   @li true, if the SystemEigen is a generalized one
   @li false otherwise
   **

   @fn SystemEigen::getEigenValuesNumber
   @return Returns the number of computed eigenvalues
   **

   @fn SystemEigen::getEigenValues
   @return Returns the computed eigenvalues
   **

   @fn SystemEigen::getEigenVectors
   @return Returns the computed eigenvectors
   **

   @fn SystemEigen::setNumberOfEigenValues
   @param nEigenValues A natural number

   Sets the number of eigenvalues computed by
   SystemEigen::solve() to the given number
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline bool SystemEigen::isGeneral(void) const{
  return general;
}

inline size_t SystemEigen::getEigenValuesNumber(void) const{
  return nEigenValues;
}

inline const std::vector<std::complex<double> >&
SystemEigen::getEigenValues(void) const{
  return *eigenValue;
}

inline const std::vector<fullVector<std::complex<double> > >&
SystemEigen::getEigenVectors(void) const{
  return *eigenVector;
}

#endif
