#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include <string>

#include "System.h"
#include "WriterMsh.h"

class Solution{
 private:
  WriterMsh* writer;

 public:
   Solution(const System& system);
  ~Solution(void);

  void write(const std::string name) const;
};

#endif
