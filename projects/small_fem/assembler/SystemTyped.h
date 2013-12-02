#ifndef _SYSTEMTYPED_H_
#define _SYSTEMTYPED_H_

#include "SystemAbstract.h"

#include "GroupOfElement.h"
#include "SolverVector.h"
#include "SolverMatrix.h"
#include "FormulationTyped.h"
#include "fullMatrix.h"

/**
   @interface SystemTyped
   @brief Type of scalar used by a linear system assembler

   This interface adds to SystemAbstract the notion of scalar type.
   A scalar may be:
   @li real
   @li complex

   This interface enables access to the finite element solution
   comuted by the SystemAbstract.
 */

template<typename scalar>
class SystemTyped: public SystemAbstract{
 protected:
  typedef
    scalar (FormulationTyped<scalar>::*formulationPtr)(size_t dofI,
                                                       size_t dofJ,
                                                       size_t elementId) const;
 protected:
  const FormulationTyped<scalar>* formulation;

 public:
  virtual ~SystemTyped(void);

  virtual std::string getType(void) const;
  virtual void getSolution(fullVector<scalar>& sol, size_t nSol) const = 0;
  virtual void getSolution(fullVector<scalar>& sol)              const = 0;

 protected:
  void assemble(SolverMatrix<scalar>& A,
                SolverVector<scalar>& b,
                size_t elementId,
                const GroupOfDof& group,
                formulationPtr& term);
};


/**
   @fn SystemTyped::~SystemTyped
   Deletes this SystemTyped
   **

   @fn SystemTyped::getSolution(fullVector<scalar>&, size_t)
   @param sol A vector
   @param nSol An integer
   Allocates and populates the given vector with the nSolth solution vector
   computed by SystemAbstract::solve()
   **

   @fn SystemTyped::getSolution(fullVector<scalar>&)
   @param sol A vector
   Same as SystemTyped::getSolution(sol, 0)
   **

   @internal
   @fn SystemTyped::assemble
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

#include "SystemTypedInclusion.h"

#endif
