#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "AbstractSystem.h"

/**
   @class System
   @brief This class assembles a linear system

   This class assembles a linear system,
   described by a Formulation
 */

class System: public AbstractSystem{
 protected:
  linearSystemPETSc<double>* linSys;
  fullVector<double>*        x;

 public:
  System(const Formulation& formulation);
  virtual ~System(void);

  const fullVector<double>& getSol(void) const;

  virtual void assemble(void);
  virtual void solve(void);
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
