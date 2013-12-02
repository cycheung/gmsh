#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "SystemTyped.h"

/**
   @class System
   @brief This class assembles a linear system

   This class assembles a linear system described by a Formulation.

   The Solver used is <a href="http://graal.ens-lyon.fr/MUMPS/index.php">MUMPS
   </a>.
 */

class System: public SystemTyped<double>{
 protected:
  SolverMatrix<double>* A;
  SolverVector<double>* b;
  fullVector<double>*   x;

 public:
  System(const FormulationTyped<double>& formulation);
  virtual ~System(void);

  virtual size_t getNComputedSolution(void)                        const;
  virtual void   getSolution(fullVector<double>& sol, size_t nSol) const;
  virtual void   getSolution(fullVector<double>& sol)              const;

  virtual void assemble(void);
  virtual void solve(void);

  virtual void   addSolution(FEMSolution& feSol) const;
  virtual void   writeMatrix(std::string fileName,
                             std::string matrixName) const;
};


/**
   @fn System::System(const Formulation&)
   @param formulation A Formulation that gives the way to assemble the system

   Instantiates a new System
   **

   @fn System::~System

   Deletes this System
*/

#endif
