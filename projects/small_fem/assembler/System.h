#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "SystemAbstract.h"

/**
   @class System
   @brief This class assembles a linear system

   This class assembles a linear system described by a Formulation.

   The Solver used is <a href="http://graal.ens-lyon.fr/MUMPS/index.php">MUMPS
   </a>.
 */
template<typename scalar>
class System: public SystemAbstract<scalar>{
 protected:
  SolverMatrix<scalar>* A;
  SolverVector<scalar>* b;
  fullVector<scalar>*   x;

 public:
  System(const Formulation<scalar>& formulation);
  virtual ~System(void);

  virtual size_t getNComputedSolution(void)                           const;
  virtual void   getSolution(fullVector<scalar>& sol, size_t nSol)    const;
  virtual void   getSolution(fullVector<scalar>& sol)                 const;
  virtual void   getSolution(std::map<Dof, scalar>& sol, size_t nSol) const;
  virtual void   getSolution(std::map<Dof, scalar>& sol)              const;
  virtual void   getSolution(FEMSolution<scalar>& feSol)              const;

  virtual void assemble(void);
  virtual void solve(void);

  virtual void writeMatrix(std::string fileName,
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

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SystemInclusion.h"


#endif
