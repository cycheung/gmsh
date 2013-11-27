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

class System: public SystemAbstract{
 protected:
  SolverMatrix*       A;
  SolverVector*       b;
  fullVector<double>* x;

 public:
  System(const Formulation& formulation);
  virtual ~System(void);

  const fullVector<double>& getSol(void) const;

  virtual void assemble(void);
  virtual void solve(void);

  virtual void addSolution(FEMSolution& feSol) const;
  virtual void writeMatrix(std::string fileName, std::string matrixName) const;

 protected:
  System(void);
};


/**
   @fn System::System(const Formulation&)
   @param formulation A Formulation that gives the way to assemble the system

   Instantiates a new System
   **

   @fn System::~System

   Deletes this System
   **

   @fn System::getSol
   @return Returns the solution of the linear system
   **
*/

/////////////////////
// Inline Function //
/////////////////////

inline const fullVector<double>& System::getSol(void) const{
  return *x;
}

#endif
