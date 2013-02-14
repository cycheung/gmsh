#ifndef _SYSTEMEIGEN_H_
#define _SYSTEMEIGEN_H_

#include "SystemAbstract.h"
#include "EigenSolver.h"

/**
   @class SystemEigen
   @brief This class assembles an Eigenvalue System

   This class assembles an Eigenvalue system,
   described by a Formulation.@n

   The Eigenvalue Problem can be @em generalized or not:
   @li An Eigenvalue Problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$
   @li An Generalized Eigenvalue Problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$

   @see Formulation::isGeneral()
 */

class SystemEigen: public SystemAbstract{
 private:
  bool general;

  linearSystemPETSc<double>* linSysA;
  linearSystemPETSc<double>* linSysB;
  EigenSolver*               eSys;

  std::vector<std::complex<double> >* eigenValue;
  std::vector<std::vector<std::complex<double> > >* eigenVector;
  unsigned int nEigenValues;

 public:
  SystemEigen(const Formulation& formulation);
  virtual ~SystemEigen(void);

  bool isGeneral(void) const;

  unsigned int                                       getEigenValuesNumber(void) const;
  const std::vector<std::complex<double> >&               getEigenValues(void)  const;
  const std::vector<std::vector<std::complex<double> > >& getEigenVectors(void) const;

  void setNumberOfEigenValues(unsigned int nEigenValues);

  virtual void assemble(void);
  virtual void solve(void);
};


/**
   @fn SystemEigen::SystemEigen
   @param formulation An Formulation that
   gives the way to assemble the Eigenvalue System

   Instantiated a new SystemEigen
   ***

   @fn SystemEigen::~SystemEigen
   Deletes this SystemEigen
   **

   @fn SystemEigen::isGeneral
   @return Returns:
   @li @c true, if the SystemEigen is a Generalized one
   @li @c false otherwise
   **

   @fn SystemEigen::getEigenValueNumber
   @return Returns the number of @em computed
   Eigenvalues
   **

   @fn SystemEigen::getEigenValues
   @return Returns the @em computed Eigenvalues
   **

   @fn SystemEigen::getEigenVectors
   @return Returns the @em computed Eigenvectors
   **

   @fn SystemEigen::setNumberOfEigenValues
   @param nEigenValues A natural number
   Set the number of eigenvalues computed by
   the solving method to the given number
*/

//////////////////////
// Inline Functions //
//////////////////////

inline bool SystemEigen::isGeneral(void) const{
  return general;
}

inline unsigned int SystemEigen::getEigenValuesNumber(void) const{
  return nEigenValues;
}

inline const std::vector<std::complex<double> >&
SystemEigen::getEigenValues(void) const{
  return *eigenValue;
}

inline const std::vector<std::vector<std::complex<double> > >&
SystemEigen::getEigenVectors(void) const{
  return *eigenVector;
}

#endif
