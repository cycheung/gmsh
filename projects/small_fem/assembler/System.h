#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "fullMatrix.h"
#include "GroupOfElement.h"
#include "GroupOfDof.h"

#include "DofManager.h"
#include "FunctionSpace.h"
#include "Formulation.h"

#include <string>

/**
   @class System
   @brief This class assembles a linear system

   This class assembles the linear system, 
   described by a Formulation 
  
   @warning
   We can @em only assemble Dof related to a MElement@n

   @todo
   Assembly of @em NON Element related Dof
 */

class System{
 private:
  bool isAssembled;

  fullMatrix<double>* A;
  fullVector<double>* b;
  fullVector<double>* x;

  int size;

  const Formulation*   formulation;
  const FunctionSpace* fs;
  DofManager*          dofM;

 public:
   System(const Formulation& formulation);
  ~System(void);

  fullMatrix<double>& getMatrix(void) const;
  fullVector<double>& getRHS(void) const;
  fullVector<double>& getSol(void) const;

  unsigned int getSize(void) const;

  const FunctionSpace& getFunctionSpace(void) const;
  const DofManager&    getDofManager(void) const;

  void fixBC(const GroupOfElement& goe, double value);
  void assemble(void);
  void solve(void);

 private:
  void assemble(GroupOfDof& group);
};


/**
   @fn System::System
   @param formulation A Formulation, 
   giving the way to assemble the system
   
   Instantiated a new System
   ***

   @fn System::~System
   Deletes this System
   **

   @fn System::getMatrix
   @return Returns the matrix of the the linear system
   **

   @fn System::getRHS
   @return Returns the Right Hand Side of the the linear system
   **

   @fn System::getSol
   @return Returns the solution of the the linear system
   **

   @fn System::getFunctionSpace
   @return Returns the FunctionSpace used by the System
   **

   @fn System::getDofManager
   @return Returns the DofManager used by the System
   **

   @fn System::fixBC
   @param goe The GroupOfElement on which the bondary condtion shall be imposed
   @param value The value of the bondary condition
   
   Fix the given Boundary Condition
   **

   @fn System::assemble
   Assembles the linear system
   **

   @fn System::solve
   Solves the linear system
   @note If the System is @em not @em assembled,@n
   the assembly method will be called
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline fullMatrix<double>& System::getMatrix(void) const{
  return *A;
}

inline fullVector<double>& System::getRHS(void) const{
  return *b;
}

inline fullVector<double>& System::getSol(void) const{
  return *x;
}

inline unsigned int System::getSize(void) const{
  return A->size1();
}

inline const FunctionSpace& System::getFunctionSpace(void) const{
  return *fs;
}

inline const DofManager& System::getDofManager(void) const{
  return *dofM;
}

#endif
