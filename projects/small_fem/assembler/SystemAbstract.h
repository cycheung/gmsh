#ifndef _SYSTEMABSTRACT_H_
#define _SYSTEMABSTRACT_H_

#include <petscmat.h>
#include <petscvec.h>

#include <vector>
#include <set>

#include "Dof.h"
#include "DofManager.h"
#include "FunctionSpace.h"
#include "Formulation.h"
#include "GroupOfElement.h"

#include "fullMatrix.h"
#include "SparseMatrix.h"

/**
   @interface SystemAbstract
   @brief Common Interface for Linear Systems Assemblers

   This is a common interface for Linear Systems Assemblers
 */

class SystemAbstract{
 protected:
  typedef
    double (Formulation::*formulationPtr)(size_t dofI,
                                          size_t dofJ,
                                          size_t elementId) const;

  typedef std::vector<std::set<size_t> > UniqueSparsity;

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

  virtual void assemble(void) = 0;
  virtual void solve(void)    = 0;

 protected:
  void assemble(SparseMatrix& A,
                fullVector<double>& b,
                size_t elementId,
                const GroupOfDof& group,
                formulationPtr& term);

  void assemble(Mat& A,
                Vec& b,
                size_t elementId,
                const GroupOfDof& group,
                formulationPtr& term);

  void assemble(Mat& A,
                size_t elementId,
                const GroupOfDof& group,
                formulationPtr& term);

  void sparsity(PetscInt* nonZero,
                UniqueSparsity& uniqueSparsity,
                const GroupOfDof& group);
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
   @return Returns the number of unknowns in the the linear system
   **

   @fn SystemAbstract::getDofManager
   @return Returns the DofManager used by the system
   **

   @fn SystemAbstract::getFunctionSpace
   @return Returns the FunctionSpace used by the system
   **

   @fn SystemAbstract::fixCoef
   @param goe A GroupOfElement
   @param value A real value

   Fixes the Coefficients (Dof%s), associated to the given
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

   @fn SystemAbstract::assemble
   Assembles the linear system
   **

   @fn SystemAbstract::solve
   Solves the linear system
   **
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
