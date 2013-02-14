#ifndef _SYSTEMINSTRUMENTED_H_
#define _SYSTEMINSTRUMENTED_H_

#include "System.h"

/**
   @class SystemInstrumented
   @brief This class is a System with Timer%s

   This class is a System with Timer%s
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

 private:
  void assemble(const GroupOfDof& group);
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
