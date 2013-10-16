#ifndef _SYSTEM_H_
#define _SYSTEM_H_

// #include <petscmat.h>
// #include <petscvec.h>

#include "SystemAbstract.h"
#include "ThreadVector.h"
#include "SparseMatrix.h"

/**
   @class System
   @brief This class assembles a linear system

   This class assembles a linear system,
   described by a Formulation
 */

class System: public SystemAbstract{
 protected:
  SparseMatrix* A;
  ThreadVector* b;

  /*
  Mat* A;
  Vec* b;
  Vec* xPetsc;
  */

  fullVector<double>* x;

 public:
  System(const Formulation& formulation);
  virtual ~System(void);

  const fullVector<double>& getSol(void) const;

  virtual void assemble(void);
  virtual void solve(void);

 protected:
  System(void);
};


/**
   @fn System::System(const Formulation&)
   @param formulation A Formulation that
   gives the way to assemble the system

   Instantiated a new System
   ***

   @fn System::~System
   Deletes this System
   **

   @fn System::getSol
   @return Returns the solution of the linear system
*/

/////////////////////
// Inline Function //
/////////////////////

inline const fullVector<double>& System::getSol(void) const{
  return *x;
}

#endif
