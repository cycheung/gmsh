#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "AbstractSystem.h"
#include "Formulation.h"

#include "linearSystemPETSc.h"

/**
   @class System
   @brief This class assembles a linear system

   This class assembles a linear system,
   described by a Formulation

   @warning
   We can @em only assemble Dof related to a MElement@n

   @todo
   Assembly of @em NON Element related Dof@n
   Allow multiple basis for dirichelt
 */

class System: public AbstractSystem{
 protected:
  const Formulation*         formulation;
  linearSystemPETSc<double>* linSys;
  fullVector<double>*        x;

 public:
  System(const Formulation& formulation);
  virtual ~System(void);

  const fullVector<double>& getSol(void) const;

  virtual void assemble(void);
  virtual void solve(void);

 private:
  void assemble(GroupOfDof& group);
  void sparsity(GroupOfDof& group);
};


/**
   @fn System::System
   @param formulation A Formulation that
   gives the way to assemble the system

   Instantiated a new System
   ***

   @fn System::~System
   Deletes this System
   **

   @fn System::getSol
   @return Returns the solution of the the linear system
*/

/////////////////////
// Inline Function //
/////////////////////

inline const fullVector<double>& System::getSol(void) const{
  return *x;
}

#endif
