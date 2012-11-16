#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "fullMatrix.h"
#include "GroupOfElement.h"
#include "GroupOfDof.h"

#include "DofManager.h"
#include "FunctionSpace.h"
#include "Formulation.h"

#include "linearSystemPETSc.h"

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
  bool assembled;
  bool solved;

  linearSystemPETSc<double>* linSys;
  fullVector<double>*        x;
  int size;

  const Formulation*   formulation;
  const FunctionSpace* fs;
  DofManager*          dofM;

 public:
   System(const Formulation& formulation);
  ~System(void);

  unsigned int              getSize(void) const;
  const fullVector<double>& getSol(void) const;

  const FunctionSpace& getFunctionSpace(void) const;
  const DofManager&    getDofManager(void) const;

  bool isAssembled(void) const;
  bool isSolved(void) const;

  void fixCoef(const GroupOfElement& goe, 
	       double value);
  
  void dirichlet(const GroupOfElement& goe, 
		 double (*f)(fullVector<double>& xyz));
  
  void dirichlet(const GroupOfElement& goe, 
		 fullVector<double> (*f)(fullVector<double>& xyz));



  void assemble(void);
  void solve(void);

 private:
  void assemble(GroupOfDof& group);
  void sparcity(GroupOfDof& group);
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

   @fn System::getSol
   @return Returns the solution of the the linear system
   **

   @fn System::getFunctionSpace
   @return Returns the FunctionSpace used by the System
   **

   @fn System::getDofManager
   @return Returns the DofManager used by the System
   **

   @fn System::fixCoef(const GroupOfElement& goe, double value)
   @param goe A GroupOfElement 
   @param value A real value
   
   Fixes the Coefficients (Dof%s) associated the the given 
   GroupOfElement to the given value

   @note These Coefficients are (Dof%s) weights of the 
   Basis Function defined on the given GroupOfElement
   **

   @fn System::fixCoef(const std::vector<std::pair<const Dof*, double> >& value);
   @param value A vector of pair of Dof and Real Value
   
   Fixes the @f$i@f$th Dof (Coefficient) to the 
   @f$i@f$th value 
   (for @f$i@f$ rangin from 0 to value.size() - 1)

   @note These Dof%s (Coefficients) are weights of the 
   Basis Function defined on the given GroupOfElement
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

inline const fullVector<double>& System::getSol(void) const{
  return *x;
}

inline unsigned int System::getSize(void) const{
  return size;
}

inline const FunctionSpace& System::getFunctionSpace(void) const{
  return *fs;
}

inline const DofManager& System::getDofManager(void) const{
  return *dofM;
}

inline bool System::isAssembled(void) const{
  return assembled;
}

inline bool System::isSolved(void) const{
  return solved;
}

#endif
