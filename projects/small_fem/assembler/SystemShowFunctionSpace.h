#ifndef _SYSTEMSHOWFUNCTIONSPACE_H_
#define _SYSTEMSHOWFUNCTIONSPACE_H_

#include "System.h"

/**
   @class SystemShowFunctionSpace
   @brief This class shows a FunctionSpace

   This class is a System shows a FunctionSpace,
   by returning the appropriate solution vector.@n

   For example, SystemShowFunctionSpace(FunctionSpace&, functionNumber),
   returns the solution vector such that: x(i)
   @li is equal to 1 if i = functionNumber
   @li is equal to 0 otherwise

   This class is a particular System.
 */

class SystemShowFunctionSpace: public System{
 protected:
  unsigned int fNumber;

 public:
  SystemShowFunctionSpace(const FunctionSpace& fs,
                          unsigned int functionNumber);

  virtual ~SystemShowFunctionSpace(void);

  virtual void assemble(void);
  virtual void solve(void);
};


/**
   @fn SystemShowFunctionSpace::SystemShowFunctionSpace
   @param fs the FunctionSpace to show
   @param functionNumber the number of the function
   (of the FunctionSpace) to show

   Instantiated a new SystemShowFunctionSpace
   ***

   @fn SystemShowFunctionSpace::~SystemShowFunctionSpace
   Deletes this SystemShowFunctionSpace
*/

#endif
