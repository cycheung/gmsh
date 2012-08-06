#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "fullMatrix.h"
#include "GroupOfElement.h"
#include "GroupOfDof.h"
#include "Writer.h"

#include "DofManager.h"
#include "FunctionSpace.h"
#include "Formulation.h"

#include <string>

/**
   @class System
   @brief This class assembles a linear system

   This class assembles the linear system, that comes from a Formulation 
   and a list of Element%s.@n
  
   @warning
   We can @em only assemble Dof related to an Element@n
 */

class System{
 private:
  bool isAssembled;

  fullMatrix<double>* A;
  fullVector<double>* b;
  fullVector<double>* x;

  int size;

  const Formulation*   formulation;
  FunctionSpace*       fs;
  DofManager*          dofM;

 public:
   System(const Formulation& formulation);
  ~System(void);

  fullMatrix<double>& getMatrix(void) const;
  fullVector<double>& getRHS(void) const;
  fullVector<double>& getSol(void) const;

  FunctionSpace&     getFunctionSpace(void) const;
  DofManager&        getDofManager(void) const;
  const Formulation& getFormulation(void) const;

  void fixBC(const GroupOfElement& goe, double value);
  void assemble(void);
  void solve(void);

 private:
  void assemble(GroupOfDof& group);
};


/**
   @fn System::System(const std::vector<MElement*>& elements, 
   const Formulation& formulation)
   @param elements A list of Element%s, giving the geomtry of the problem to solve
   @param formulation A Formulation, giving the way to assemble the system
   @return A new System
 
   @fn System::~System(void)
   @return Deletes the System

   @fn Matrix& System::getMatrix(void) const
   @return Returns the assembled matrix of the the linear system

   @fn fullVector<double>& System::getRHS(void) const
   @return Returns the assembled Right Hand Side of the the linear system

   @fn fullVector<double>& System::getSol(void) const
   @return Returns the solution of the the linear system

   @fn void System::fixBC(const int physicalId, const double value)
   @param physicalId The physical @c ID on which the bondary condtion shall be imposed
   @param value The value of the bondary condition
   @return Fix a Boundary Condition on the linear system

   @fn void System::assemble(void)
   @return Assembles the linear system

   @fn void System::solve(void)
   @return Solves the linear system
   @note If the System is @em not @em assembled,@n
   the assembly method will be called
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

inline FunctionSpace& System::getFunctionSpace(void) const{
  return *fs;
}

inline DofManager& System::getDofManager(void) const{
  return *dofM;
}

inline const Formulation& System::getFormulation(void) const{
  return *formulation;
}

#endif
