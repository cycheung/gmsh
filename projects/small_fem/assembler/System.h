#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "fullMatrix.h"
#include "GroupOfElement.h"
#include "GroupOfDof.h"

#include "DofManager.h"
#include "FunctionSpace.h"
#include "Formulation.h"

#include "linearSystemPETSc.h"

#include <set>
#include <string>

/**
   @class System
   @brief This class assembles a linear system

   This class assembles a linear system,
   described by a Formulation

   @warning
   We can @em only assemble Dof related to a MElement@n

   @todo
   Assembly of @em NON Element related Dof@n
   Allow multiple basis for dirichelt
 */

class System{
 private:
  bool assembled;
  bool solved;

  std::set<const Dof*, DofComparator>* fixedOnes;

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

  const Formulation& getFormulation(void) const;
  const DofManager&  getDofManager(void) const;

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
   @param formulation A Formulation that
   gives the way to assemble the system

   Instantiated a new System
   ***

   @fn System::~System
   Deletes this System
   **

   @fn System::getSize
   @return Returns the number of @em unknowns
   of the the linear system
   **

   @fn System::getSol
   @return Returns the solution of the the linear system
   **

   @fn System::getFormulation
   @return Returns the Formulation used by the System
   **

   @fn System::getDofManager
   @return Returns the DofManager used by the System
   **

   @fn System::isAssembled
   @return Returns:
   @li @c true, if the System has been assembled
   @li @c false otherwise
   **

   @fn System::isSolved
   @return Returns:
   @li @c true, if the System has been solved
   @li @c false otherwise
   **

   @fn System::fixCoef
   @param goe A GroupOfElement
   @param value A real value

   Fixes the Coefficients (Dof%s), associated the the given
   GroupOfElement, to the given value

   @note These Coefficients are (Dof%s) weights of the
   Basis Function defined on the given GroupOfElement
   **

   @fn System::dirichlet(const GroupOfElement& goe,  double (*f)(fullVector<double>& xyz));
   @param goe A GroupOfElement
   @param f A scalar Function

   Imposes a @em scalar Dirichlet Condition (given by @c f) on the
   given GroupOfElement
   **

   @fn System::dirichlet(const GroupOfElement& goe, fullVector<double> (*f)(fullVector<double>& xyz));
   @param goe A GroupOfElement
   @param f A vectorial Function

   Imposes a @em vectorial Dirichlet Condition (given by @c f) on the
   given GroupOfElement
   **

   @fn System::assemble
   Assembles the linear system
   **

   @fn System::solve
   Solves the linear system
   @note If the System is @em not @em assembled,
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

inline const Formulation& System::getFormulation(void) const{
  return *formulation;
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
