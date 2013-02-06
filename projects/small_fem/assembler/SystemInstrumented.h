#ifndef _SYSTEMINSTRUMENTED_H_
#define _SYSTEMINSTRUMENTED_H_

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
   @class SystemInstrumented
   @brief This class assembles a linear system (with Timer%s)

   This class assembles a linear system,
   described by a Formulation.

   This class got also Timer%s.

   @warning
   We can @em only assemble Dof related to a MElement@n

   @todo
   Assembly of @em NON Element related Dof@n
   Allow multiple basis for dirichelt
 */

class SystemInstrumented{
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
  double preAlloc;

  double dofLookTime;
  int    dofLookCall;

  double totLHSTime;
  int    totLHSCall;

  double totAddLHSTime;
  int    totAddLHSCall;

  double totRHSTime;
  int    totRHSCall;

  double totAddRHSTime;
  int    totAddRHSCall;

 public:
   SystemInstrumented(const Formulation& formulation);
  ~SystemInstrumented(void);

  unsigned int              getSize(void) const;
  const fullVector<double>& getSol(void) const;

  const Formulation& getFormulation(void) const;
  const DofManager&  getDofManager(void) const;

  bool isAssembled(void) const;
  bool isSolved(void) const;

  void fixCoef(const GroupOfElement& goe,
	       double value);

  void dirichlet(GroupOfElement& goe,
		 double (*f)(fullVector<double>& xyz));

  void dirichlet(GroupOfElement& goe,
		 fullVector<double> (*f)(fullVector<double>& xyz));

  void assemble(void);
  void solve(void);

 private:
  void assemble(GroupOfDof& group);
  void sparsity(GroupOfDof& group);
};


/**
   @fn SystemInstrumented::SystemInstrumented
   @param formulation A Formulation that
   gives the way to assemble the system

   Instantiated a new SystemInstrumented
   ***

   @fn SystemInstrumented::~SystemInstrumented
   Deletes this SystemInstrumented
   **

   @fn SystemInstrumented::getSize
   @return Returns the number of @em unknowns
   of the the linear system
   **

   @fn SystemInstrumented::getSol
   @return Returns the solution of the the linear system
   **

   @fn SystemInstrumented::getFormulation
   @return Returns the Formulation used by the System
   **

   @fn SystemInstrumented::getDofManager
   @return Returns the DofManager used by the System
   **

   @fn SystemInstrumented::isAssembled
   @return Returns:
   @li @c true, if the System has been assembled
   @li @c false otherwise
   **

   @fn SystemInstrumented::isSolved
   @return Returns:
   @li @c true, if the System has been solved
   @li @c false otherwise
   **

   @fn SystemInstrumented::fixCoef
   @param goe A GroupOfElement
   @param value A real value

   Fixes the Coefficients (Dof%s), associated the the given
   GroupOfElement, to the given value

   @note These Coefficients are (Dof%s) weights of the
   Basis Function defined on the given GroupOfElement
   **

   @fn SystemInstrumented::dirichlet(const GroupOfElement& goe,  double (*f)(fullVector<double>& xyz));
   @param goe A GroupOfElement
   @param f A scalar Function

   Imposes a @em scalar Dirichlet Condition (given by @c f) on the
   given GroupOfElement
   **

   @fn SystemInstrumented::dirichlet(const GroupOfElement& goe, fullVector<double> (*f)(fullVector<double>& xyz));
   @param goe A GroupOfElement
   @param f A vectorial Function

   Imposes a @em vectorial Dirichlet Condition (given by @c f) on the
   given GroupOfElement
   **

   @fn SystemInstrumented::assemble
   Assembles the linear system
   **

   @fn SystemInstrumented::solve
   Solves the linear system
   @note If the System is @em not @em assembled,
   the assembly method will be called
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline const fullVector<double>& SystemInstrumented::getSol(void) const{
  return *x;
}

inline unsigned int SystemInstrumented::getSize(void) const{
  return size;
}

inline const Formulation& SystemInstrumented::getFormulation(void) const{
  return *formulation;
}

inline const DofManager& SystemInstrumented::getDofManager(void) const{
  return *dofM;
}

inline bool SystemInstrumented::isAssembled(void) const{
  return assembled;
}

inline bool SystemInstrumented::isSolved(void) const{
  return solved;
}

#endif
