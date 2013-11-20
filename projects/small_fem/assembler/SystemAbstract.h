#ifndef _SYSTEMABSTRACT_H_
#define _SYSTEMABSTRACT_H_

#include <string>

#include "DofManager.h"
#include "FunctionSpace.h"
#include "Formulation.h"
#include "GroupOfElement.h"

#include "fullMatrix.h"
#include "SolverVector.h"
#include "SolverMatrix.h"

/**
   @interface SystemAbstract
   @brief Common interface for linear systems assemblers

   This is a common interface for linear systems assemblers
 */

class SystemAbstract{
 protected:
  typedef double (Formulation::*formulationPtr)(size_t dofI,
                                                size_t dofJ,
                                                size_t elementId) const;
 protected:
  bool assembled;
  bool solved;

  const Formulation*   formulation;
  const FunctionSpace* fs;
  DofManager*          dofM;

 public:
  virtual ~SystemAbstract(void);

  bool isAssembled(void) const;
  bool isSolved(void)    const;

  size_t               getSize(void)          const;
  const DofManager&    getDofManager(void)    const;
  const FunctionSpace& getFunctionSpace(void) const;

  void fixCoef(const GroupOfElement& goe,
               double value);

  void dirichlet(GroupOfElement& goe,
                 double (*f)(fullVector<double>& xyz));

  void dirichlet(GroupOfElement& goe,
                 fullVector<double> (*f)(fullVector<double>& xyz));

  virtual void writeMatrix(std::string fileName, std::string matrixName) const;

  virtual void assemble(void) = 0;
  virtual void solve(void)    = 0;

 protected:
  void assemble(SolverMatrix& A,
                SolverVector& b,
                size_t elementId,
                const GroupOfDof& group,
                formulationPtr& term);
};


/**
   @fn SystemAbstract::~SystemAbstract
   Deletes this SystemAbstract
   **

   @fn SystemAbstract::isAssembled
   @return Returns:
   @li true, if the system has been assembled
   @li false otherwise
   **

   @fn SystemAbstract::isSolved
   @return Returns:
   @li true, if the system has been solved
   @li false otherwise
   **

   @fn SystemAbstract::getSize
   @return Returns the number of unknowns in this linear system
   **

   @fn SystemAbstract::getDofManager
   @return Returns the DofManager used by this system
   **

   @fn SystemAbstract::getFunctionSpace
   @return Returns the FunctionSpace used by this system
   **

   @fn SystemAbstract::fixCoef
   @param goe A GroupOfElement
   @param value A real value

   Fixes the coefficients (Dof%s), associated to the given
   GroupOfElement, to the given value
   **

   @fn SystemAbstract::dirichlet(GroupOfElement& goe,  double (*f)(fullVector<double>& xyz));
   @param goe A GroupOfElement
   @param f A scalar Function

   Imposes a scalar Dirichlet Condition (given by f) on the given GroupOfElement
   **

   @fn SystemAbstract::dirichlet(GroupOfElement& goe, fullVector<double> (*f)(fullVector<double>& xyz));
   @param goe A GroupOfElement
   @param f A vectorial Function

   Imposes a vectorial Dirichlet Condition (given by @c f) on the
   given GroupOfElement
   **

   @fn SystemAbstract::writeMatrix
   @param fileName A string
   @param matrixName A string

   Writes this system matrix in Octave/Matlab format, with the given name,
   into the given file
   **

   @fn SystemAbstract::assemble(void)
   Assembles this linear system
   **

   @fn SystemAbstract::solve(void)
   Solves this linear system
   **

   @internal
   @fn SystemAbstract::assemble(SolverMatrix&, SolverVector&, size_t, const GroupOfDof&, formulationPtr&)
   @param A The matrix to assemble
   @param b The right hand side to assemble
   @param elementId The mesh Element ID to assemble
   @param group The GroupOfDof to assemble
   @param term The Formulation to use in the assembly

   Assembles the given values
   @endinternal
*/

//////////////////////
// Inline Functions //
//////////////////////

inline bool SystemAbstract::isAssembled(void) const{
  return assembled;
}

inline bool SystemAbstract::isSolved(void) const{
  return solved;
}

inline size_t SystemAbstract::getSize(void) const{
  return dofM->getUnfixedDofNumber();
}

inline const DofManager& SystemAbstract::
getDofManager(void) const{
  return *dofM;
}

inline const FunctionSpace& SystemAbstract::
getFunctionSpace(void) const{
  return *fs;
}

#endif
