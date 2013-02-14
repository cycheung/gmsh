#ifndef _ABSTRACTSYSTEM_H_
#define _ABSTRACTSYSTEM_H_

#include "Dof.h"
#include "DofManager.h"
#include "FunctionSpace.h"
#include "Formulation.h"
#include "GroupOfElement.h"
#include "fullMatrix.h"

#include "linearSystemPETSc.h"

/**
   @interface Common Interface for Linear Systems
   @brief Common Interface for Linear Systems Assemblers

   This is a common Interface for Linear Systems Assemblers

   @warning
   We can @em only assemble Dof related to a MElement@n

   @todo
   Assembly of @em NON Element related Dof@n
   Allow multiple basis for dirichelt
 */

class AbstractSystem{
 protected:
  typedef
    double (Formulation::*formulationPtr)(unsigned int dofI,
                                          unsigned int dofJ,
                                          const GroupOfDof& god) const;

 protected:
  bool assembled;
  bool solved;

  const Formulation*   formulation;
  const FunctionSpace* fs;
  DofManager*          dofM;

 public:
  virtual ~AbstractSystem(void);

  bool isAssembled(void) const;
  bool isSolved(void)    const;

  unsigned int         getSize(void)          const;
  const DofManager&    getDofManager(void)    const;
  const FunctionSpace& getFunctionSpace(void) const;

  void fixCoef(const GroupOfElement& goe,
               double value);

  void dirichlet(GroupOfElement& goe,
                 double (*f)(fullVector<double>& xyz));

  void dirichlet(GroupOfElement& goe,
                 fullVector<double> (*f)(fullVector<double>& xyz));

  virtual void assemble(void) = 0;
  virtual void solve(void)    = 0;

 protected:
  void assemble(linearSystemPETSc<double>& sys,
                GroupOfDof& group,
                formulationPtr& term);

  void sparsity(linearSystemPETSc<double>& sys,
                GroupOfDof& group);
};


/**
   @fn AbstractSystem::~AbstractSystem
   Deletes this AbstractSystem
   **

   @fn AbstractSystem::isAssembled
   @return Returns:
   @li @c true, if the system has been assembled
   @li @c false otherwise
   **

   @fn AbstractSystem::isSolved
   @return Returns:
   @li @c true, if the system has been solved
   @li @c false otherwise
   **

   @fn AbstractSystem::getSize
   @return Returns the number of @em unknowns
   of the the linear system
   **

   @fn AbstractSystem::getDofManager
   @return Returns the DofManager used by the system
   **

   @fn AbstractSystem::getFunctionSpace
   @return Returns the FunctionSpace used by the system
   **

   @fn AbstractSystem::fixCoef
   @param goe A GroupOfElement
   @param value A real value

   Fixes the Coefficients (Dof%s), associated the the given
   GroupOfElement, to the given value

   @note These Coefficients are (Dof%s) weights of the
   Basis Function defined on the given GroupOfElement
   **

   @fn AbstractSystem::dirichlet(const GroupOfElement& goe,  double (*f)(fullVector<double>& xyz));
   @param goe A GroupOfElement
   @param f A scalar Function

   Imposes a @em scalar Dirichlet Condition (given by @c f) on the
   given GroupOfElement
   **

   @fn AbstractSystem::dirichlet(const GroupOfElement& goe, fullVector<double> (*f)(fullVector<double>& xyz));
   @param goe A GroupOfElement
   @param f A vectorial Function

   Imposes a @em vectorial Dirichlet Condition (given by @c f) on the
   given GroupOfElement
   **

   @fn AbstractSystem::assemble
   Assembles the linear system
   **

   @fn AbstractSystem::solve
   Solves the linear system
   @note If the AbstractSystem is @em not @em assembled,
   the assembly method will be called
*/

//////////////////////
// Inline Functions //
//////////////////////

inline bool AbstractSystem::isAssembled(void) const{
  return assembled;
}

inline bool AbstractSystem::isSolved(void) const{
  return solved;
}

inline unsigned int AbstractSystem::getSize(void) const{
  return fs->dofNumber();
}

inline const DofManager& AbstractSystem::
getDofManager(void) const{
  return *dofM;
}

inline const FunctionSpace& AbstractSystem::
getFunctionSpace(void) const{
  return *fs;
}

#endif
