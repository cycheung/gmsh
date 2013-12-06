#ifndef _SYSTEMABSTRACT_H_
#define _SYSTEMABSTRACT_H_

#include <string>

#include "DofManager.h"
#include "Formulation.h"
#include "FunctionSpace.h"
#include "GroupOfElement.h"
#include "FEMSolution.h"

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "fullMatrix.h"

/**
   @interface SystemAbstract
   @brief Common interface for linear systems assemblers

   This is a common interface for linear systems assemblers.
   A SystemAbstract may be of multiple scalar type.
 */

template<typename scalar>
class SystemAbstract{
 protected:
  typedef
    scalar (Formulation<scalar>::*formulationPtr)(size_t dofI,
                                                  size_t dofJ,
                                                  size_t elementId) const;
 protected:
  static const scalar minusSign;

 protected:
  bool assembled;
  bool solved;

  const Formulation<scalar>* formulation;
  const FunctionSpace*       fs;
  DofManager<scalar>*        dofM;

 public:
  virtual ~SystemAbstract(void);

  bool isAssembled(void) const;
  bool isSolved(void)    const;

  size_t                    getSize(void)          const;
  const DofManager<scalar>& getDofManager(void)    const;
  const FunctionSpace&      getFunctionSpace(void) const;

  void constraint(const std::map<Dof, scalar>& constr);

  virtual void assemble(void) = 0;
  virtual void solve(void)    = 0;

  virtual size_t getNComputedSolution(void)                           const = 0;
  virtual void   getSolution(fullVector<scalar>& sol, size_t nSol)    const = 0;
  virtual void   getSolution(fullVector<scalar>& sol)                 const = 0;
  virtual void   getSolution(std::map<Dof, scalar>& sol, size_t nSol) const = 0;
  virtual void   getSolution(std::map<Dof, scalar>& sol)              const = 0;
  virtual void   getSolution(FEMSolution<scalar>& feSol)              const = 0;

  virtual void  writeMatrix(std::string fileName, std::string matrixName) const;

 protected:
  void assemble(SolverMatrix<scalar>& A,
                SolverVector<scalar>& b,
                size_t elementId,
                const GroupOfDof& group,
                formulationPtr& term,
                const Formulation<scalar>& formulation);
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
   @param constr A map of Dof%s and scalar
   Constraints this SystemAbstract with the given Dof%s
   and their associated values
   **

   @fn SystemAbstract::assemble(void)
   Assembles this linear system
   **

   @fn SystemAbstract::solve(void)
   Solves this linear system
   **

   @fn SystemAbstract::getNComputedSolution(void)
   @return The number of computed solution by SystemAbstract::solve()
   **

   @fn SystemAbstract::getSolution(fullVector<scalar>&, size_t)
   @param sol A vector
   @param nSol An integer
   Allocates and populates the given vector with the nSolth solution vector
   computed by SystemAbstract::solve()
   **

   @fn SystemAbstract::getSolution(fullVector<scalar>&)
   @param sol A vector
   Same as SystemAbstract::getSolution(sol, 0)
   **

   @fn SystemAbstract::getSolution(std::map<Dof, scalar>&, size_t)
   @param sol A map mapping a Dof to a scalar
   @param nSol An integer
   Fills to given map with the given pairs
   (Dof, solution associated to this Dof)
   for the nSolth solution of this SystemAbstract.
   **

   @fn SystemAbstract::getSolution(std::map<Dof, scalar>&)
   @param sol A map mapping a Dof to a scalar
   Same as SystemAbstract::getSolution(sol, 0)
   **

   @fn SystemAbstract::getSolution(FEMSolution& feSol)
   @param feSol A FEMSolution

   Adds to the given FEMSolution the computed finite element solutions.
   If no solution has been computed, and Exception is throw.
   **

   @fn SystemAbstract::writeMatrix
   @param fileName A string
   @param matrixName A string

   Writes this system matrix in Octave/Matlab format, with the given name,
   into the given file

   @internal
   @fn SystemAbstract::assemble
   @param A The matrix to assemble
   @param b The right hand side to assemble
   @param elementId The mesh Element ID to assemble
   @param group The GroupOfDof to assemble
   @param term The Formulation to use in the assembly

   Assembles the given values
   @endinternal
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SystemAbstractInclusion.h"

#endif
