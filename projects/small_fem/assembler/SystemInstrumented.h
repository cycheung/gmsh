#ifndef _SYSTEMINSTRUMENTED_H_
#define _SYSTEMINSTRUMENTED_H_

#include "System.h"

/**
   @class SystemInstrumented
   @brief This class assembles a linear system (with Timer%s)

   This class assembles a linear system,
   described by a Formulation.

   This class got also Timer%s.

   @warning
   We can @em only assemble Dof related to a MElement@n

   @todo
   Assembly of @em NON Element related Dof@n
   Allow multiple basis for dirichelt
 */

class SystemInstrumented: public System{
 public:
  double preAlloc;

  double dofLookTime;
  int    dofLookCall;

  double totLHSTime;
  int    totLHSCall;

  double totAddLHSTime;
  int    totAddLHSCall;

  double totRHSTime;
  int    totRHSCall;

  double totAddRHSTime;
  int    totAddRHSCall;

 public:
  SystemInstrumented(const Formulation& formulation);
  virtual ~SystemInstrumented(void);

  virtual void assemble(void);
  virtual void solve(void);

 private:
  void assemble(GroupOfDof& group);
  void sparsity(GroupOfDof& group);
};


/**
   @fn SystemInstrumented::SystemInstrumented
   @param formulation A Formulation that
   gives the way to assemble the system

   Instantiated a new SystemInstrumented
   ***

   @fn SystemInstrumented::~SystemInstrumented
   Deletes this SystemInstrumented
*/

#endif
