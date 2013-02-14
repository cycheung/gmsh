#ifndef _EIGENSYSTEM_H_
#define _EIGENSYSTEM_H_

#include "AbstractSystem.h"
#include "EigenSolver.h"

/**
   @class EigenSystem
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

class EigenSystem: public AbstractSystem{
 private:
  bool general;

  linearSystemPETSc<double>* linSysA;
  linearSystemPETSc<double>* linSysB;
  EigenSolver*               eSys;

  std::vector<std::complex<double> >* eigenValue;
  std::vector<std::vector<std::complex<double> > >* eigenVector;
  unsigned int nEigenValues;

 public:
  EigenSystem(const Formulation& formulation);
  virtual ~EigenSystem(void);

  bool isGeneral(void) const;

  unsigned int                                       getEigenValuesNumber(void) const;
  const std::vector<std::complex<double> >&               getEigenValues(void)  const;
  const std::vector<std::vector<std::complex<double> > >& getEigenVectors(void) const;

  void setNumberOfEigenValues(unsigned int nEigenValues);

  virtual void assemble(void);
  virtual void solve(void);
};


/**
   @fn EigenSystem::EigenSystem
   @param formulation An Formulation that
   gives the way to assemble the Eigenvalue System

   Instantiated a new EigenSystem
   ***

   @fn EigenSystem::~EigenSystem
   Deletes this EigenSystem
   **

   @fn EigenSystem::isGeneral
   @return Returns:
   @li @c true, if the EigenSystem is a Generalized one
   @li @c false otherwise
   **

   @fn EigenSystem::getEigenValueNumber
   @return Returns the number of @em computed
   Eigenvalues
   **

   @fn EigenSystem::getEigenValues
   @return Returns the @em computed Eigenvalues
   **

   @fn EigenSystem::getEigenVectors
   @return Returns the @em computed Eigenvectors
   **

   @fn EigenSystem::setNumberOfEigenValues
   @param nEigenValues A natural number
   Set the number of eigenvalues computed by
   the solving method to the given number
*/

//////////////////////
// Inline Functions //
//////////////////////

inline bool EigenSystem::isGeneral(void) const{
  return general;
}

inline unsigned int EigenSystem::getEigenValuesNumber(void) const{
  return nEigenValues;
}

inline const std::vector<std::complex<double> >&
EigenSystem::getEigenValues(void) const{
  return *eigenValue;
}

inline const std::vector<std::vector<std::complex<double> > >&
EigenSystem::getEigenVectors(void) const{
  return *eigenVector;
}

#endif
