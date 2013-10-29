#ifndef _SYSTEMEIGEN_H_
#define _SYSTEMEIGEN_H_

#include <complex>
#include "SystemAbstract.h"

/**
   @class SystemEigen
   @brief This class assembles an Eigenvalue System

   This class assembles an Eigenvalue system,
   described by a Formulation.

   The Eigenvalue problem can be generalized or not:
   @li An Eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$
   @li A Generalized Eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$
 */

class SystemEigen: public SystemAbstract{
 private:
  bool general;
  /*
  Mat* A;
  Mat* B;

  PetscInt nEigenValues;
  */
  size_t nEigenValues;
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
   gives the way to assemble the Eigenvalue System

   Instantiated a new SystemEigen
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
   @return Returns the number of computed Eigenvalues
   **

   @fn SystemEigen::getEigenValues
   @return Returns the computed Eigenvalues
   **

   @fn SystemEigen::getEigenVectors
   @return Returns the computed Eigenvectors
   **

   @fn SystemEigen::setNumberOfEigenValues
   @param nEigenValues A natural number
   Sets the number of eigenvalues computed by
   the solving method to the given number
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
