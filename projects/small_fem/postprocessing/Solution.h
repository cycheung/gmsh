#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include "System.h"
#include "WriterMsh.h"

class Solution{
 private:

 public:
  Solution(const System& system);
  virtual ~Solution(void);
};

#endif
