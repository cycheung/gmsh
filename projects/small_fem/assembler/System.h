#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "SystemAbstract.h"

/**
   @class System
   @brief This class assembles a linear system

   This class assembles a linear system,
   described by a Formulation
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
  void writeMatrix(std::string fileName, std::string matrixName) const;

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
   **

   @fn System::writeMatrix
   @param fileName A string

   Writes this System matrix in Octave/Matlab format into the given file
*/

/////////////////////
// Inline Function //
/////////////////////

inline const fullVector<double>& System::getSol(void) const{
  return *x;
}

inline void System::writeMatrix(std::string fileName,
                                std::string matrixName) const{
  A->writeToMatlabFile(fileName, matrixName);
}

#endif
