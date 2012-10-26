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

#include <string>

/**
   @class EigenSystem
   @brief This class assembles an Eigenvalue System

   This class assembles an Eigenvalue ystem, 
   described by a EigenFormulation 
  
   @warning
   We can @em only assemble Dof related to a MElement@n

   @todo
   Assembly of @em NON Element related Dof
 */

class EigenSystem{
 private:
  bool isAssembled;

  linearSystemPETSc<double>* linSysA;
  linearSystemPETSc<double>* linSysB;
  EigenSolver*               eSys;

  int size;

  const EigenFormulation* eFormulation;
  const FunctionSpace*    fs;
  DofManager*             dofM;

 public:
   EigenSystem(const EigenFormulation& eFormulation);
  ~EigenSystem(void);

  unsigned int        getSize(void) const;
  //fullVector<double>& getSol(void) const;

  const FunctionSpace& getFunctionSpace(void) const;
  const DofManager&    getDofManager(void) const;

  void fixDof(const GroupOfElement& goe, double value);
  void assemble(void);
  void solve(void);

 private:
  void assemble(GroupOfDof& group);
  void sparcity(GroupOfDof& group);
};


/**
   @fn EigenSystem::EigenSystem
   @param formulation An EigenFormulation, 
   giving the way to assemble the Eigenvalue System
   
   Instantiated a new EigenSystem
   ***

   @fn EigenSystem::~EigenSystem
   Deletes this EigenSystem
   **

   @fn EigenSystem::getSol
   @return Returns the solution of the the Eigenvalue System
   **

   @fn EigenSystem::getFunctionSpace
   @return Returns the FunctionSpace used by the Eigenvalue System
   **

   @fn EigenSystem::getDofManager
   @return Returns the DofManager used by the Eigenvalue System
   **

   @fn EigenSystem::fixDof
   @param goe A GroupOfElement 
   @param value A real value
   
   Fixes the Dof%s associated the the given 
   GroupOfElement to the given value
   **

   @fn EigenSystem::assemble
   Assembles the Eigenvalue System
   **

   @fn EigenSystem::solve
   Solves the Eigenvalue System
   @note If the EigenSystem is @em not @em assembled,@n
   the assembly method will be called
   **
*/

//////////////////////
// Inline Functions //
//////////////////////
/*
inline fullVector<double>& EigenSystem::getSol(void) const{
  return *x;
}
*/
inline unsigned int EigenSystem::getSize(void) const{
  return size;
}

inline const FunctionSpace& EigenSystem::getFunctionSpace(void) const{
  return *fs;
}

inline const DofManager& EigenSystem::getDofManager(void) const{
  return *dofM;
}

#endif
