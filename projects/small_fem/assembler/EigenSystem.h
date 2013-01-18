#ifndef _EIGENSYSTEM_H_
#define _EIGENSYSTEM_H_

#include "fullMatrix.h"
#include "GroupOfElement.h"
#include "GroupOfDof.h"

#include "DofManager.h"
#include "FunctionSpace.h"
#include "EigenFormulation.h"

#include "linearSystemPETSc.h"
#include "EigenSolver.h"

#include <set>
#include <vector>
#include <string>

/**
   @class EigenSystem
   @brief This class assembles an Eigenvalue System

   This class assembles an Eigenvalue system,
   described by a EigenFormulation.@n

   The Eigenvalue Problem can be @em generalized or not:
   @li An Eigenvalue Problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$
   @li An Generalized Eigenvalue Problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$

   @see EigenFormulation


   @warning
   We can @em only assemble Dof related to a MElement@n

   @todo
   Assembly of @em NON Element related Dof@n
   Allow multiple basis for dirichelt
 */

class EigenSystem{
 private:
  bool assembled;
  bool solved;
  bool general;

  std::set<const Dof*, DofComparator>* fixedOnes;

  linearSystemPETSc<double>* linSysA;
  linearSystemPETSc<double>* linSysB;
  EigenSolver*               eSys;

  unsigned int size;

  std::vector<std::complex<double> >* eigenValue;
  std::vector<std::vector<std::complex<double> > >* eigenVector;
  unsigned int nEigenValue;

  const EigenFormulation* eFormulation;
  const FunctionSpace*    fs;
  DofManager*             dofM;

 public:
   EigenSystem(const EigenFormulation& eFormulation);
  ~EigenSystem(void);

  unsigned int getSize(void) const;
  unsigned int getEigenValueNumber(void) const;

  const std::vector<std::complex<double> >&               getEigenValues(void) const;
  const std::vector<std::vector<std::complex<double> > >& getEigenVectors(void) const;

  const EigenFormulation& getFormulation(void) const;
  const DofManager&       getDofManager(void) const;

  bool isAssembled(void) const;
  bool isSolved(void) const;
  bool isGeneral(void) const;

  void fixCoef(const GroupOfElement& goe, double value);
  void dirichlet(const GroupOfElement& goe,
		 double (*f)(fullVector<double>& xyz));

  void assemble(void);
  void solve(unsigned int nEigenValues);

 private:
  void assemble(GroupOfDof& group);
  void assembleGeneral(GroupOfDof& group);
  void sparcity(GroupOfDof& group);
  void sparcityGeneral(GroupOfDof& group);
};


/**
   @fn EigenSystem::EigenSystem
   @param eFormulation An EigenFormulation that
   gives the way to assemble the Eigenvalue System

   Instantiated a new EigenSystem
   ***

   @fn EigenSystem::~EigenSystem
   Deletes this EigenSystem
   **

   @fn EigenSystem::getSize
   @return Returns the number of @em unknowns
   of the the Eigenvalue System
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

   @fn EigenSystem::getFormulation
   @return Returns the EigenFormulation used by the Eigenvalue System
   **

   @fn EigenSystem::getDofManager
   @return Returns the DofManager used by the Eigenvalue System
   **

   @fn EigenSystem::fixCoef(const GroupOfElement& goe, double value)
   @param goe A GroupOfElement
   @param value A real value

   Fixes the Coefficients (Dof%s) associated the the given
   GroupOfElement to the given value

   @note These Coefficients are (Dof%s) weights of the
   Basis Function defined on the given GroupOfElement
   **

   @fn EigenSystem::isAssembled
   @return Returns:
   @li @c true, if the EigenSystem has been assembled
   @li @c false otherwise
   **

   @fn EigenSystem::isSolved
   @return Returns:
   @li @c true, if the EigenSystem has been solved
   @li @c false otherwise

   @fn EigenSystem::assemble
   Assembles the Eigenvalue System
   **

   @fn EigenSystem::isGeneral
   @return Returns:
   @li @c true, if the EigenSystem is a Generalized one
   @li @c false otherwise
   **

   @fn EigenSystem::assemble
   Assembles the Eigenvalue System
   **

   @fn EigenSystem::solve
   @param nEigenValues A natural number

   Solves the Eigenvalue System@n
   The EigenSystem will compute @c nEigenValues Eigenvalues

   @note If the EigenSystem is @em not @em assembled,@n
   the assembly method will be called
*/

//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int EigenSystem::getSize(void) const{
  return size;
}

inline unsigned int EigenSystem::getEigenValueNumber(void) const{
  return nEigenValue;
}

inline const std::vector<std::complex<double> >&
EigenSystem::getEigenValues(void) const{
  return *eigenValue;
}

inline const std::vector<std::vector<std::complex<double> > >&
EigenSystem::getEigenVectors(void) const{
  return *eigenVector;
}

inline const EigenFormulation& EigenSystem::getFormulation(void) const{
  return *eFormulation;
}

inline const DofManager& EigenSystem::getDofManager(void) const{
  return *dofM;
}

inline bool EigenSystem::isAssembled(void) const{
  return assembled;
}

inline bool EigenSystem::isSolved(void) const{
  return solved;
}

inline bool EigenSystem::isGeneral(void) const{
  return general;
}

#endif
