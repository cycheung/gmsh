#ifndef _SYSTEMABSTRACT_H_
#define _SYSTEMABSTRACT_H_

#include <string>

#include "DofManager.h"
#include "Formulation.h"
#include "FunctionSpace.h"
#include "GroupOfElement.h"
#include "FEMSolution.h"

/**
   @interface SystemAbstract
   @brief Common interface for linear systems assemblers

   This is a common interface for linear systems assemblers.
   A SystemAbstract may be of multiple scalar type:
   see SystemTyped for more informations.
 */

class SystemAbstract{
 protected:
  bool assembled;
  bool solved;

  const FunctionSpace* fs;
  DofManager*          dofM;

 public:
  virtual ~SystemAbstract(void);

  bool isAssembled(void) const;
  bool isSolved(void)    const;

  size_t               getSize(void)          const;
  const DofManager&    getDofManager(void)    const;
  const FunctionSpace& getFunctionSpace(void) const;

  void constraint(const Formulation& formulation);

  virtual void assemble(void) = 0;
  virtual void solve(void)    = 0;

  virtual std::string getType(void)                       const = 0;
  virtual size_t      getNComputedSolution(void)          const = 0;
  virtual void        addSolution(FEMSolution& feSol)     const = 0;
  virtual void        writeMatrix(std::string fileName,
                                  std::string matrixName) const;
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

   @fn SystemAbstract::constraint
   @param formulation A Formulation
   Assembles and solves the given Formulation.
   The resulting solution will be used as a constraint for this AbstractSystem.
   **

   @fn SystemAbstract::assemble(void)
   Assembles this linear system
   **

   @fn SystemAbstract::solve(void)
   Solves this linear system
   **

   @fn SystemAbstract::getType
   @return Returns a string specifying of which type this SystemAbstract is
   @see SystemTyped
   **

   @fn SystemAbstract::getNComputedSolution(void)
   @return The number of computed solution by SystemAbstract::solve()

   **
   @fn SystemAbstract::addSolution(FEMSolution& feSol)
   @param feSol A FEMSolution

   Adds to the given FEMSolution the computed finite element solutions.
   If no solution has been computed, and Exception is throw.
   **

   @fn SystemAbstract::writeMatrix
   @param fileName A string
   @param matrixName A string

   Writes this system matrix in Octave/Matlab format, with the given name,
   into the given file
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
