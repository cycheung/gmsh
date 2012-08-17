#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include <string>
#include <vector>

#include "System.h"
#include "Writer.h"
#include "fullMatrix.h"

#include "Mesh.h"
#include "DofManager.h"
#include "FunctionSpace.h"

/**
   @class Solution
   @brief Writes the solution of a System

   This class can @em write the solution of a System
   into a file.@n

   The file format is defined by a Writer.
 */

class Solution{
 private:
  fullVector<double>* sol;

  const Mesh*          mesh;
  const DofManager*    dofM;
  const FunctionSpace* fs;

  int fsType;
  
  bool                              isScalar;
  std::vector<double>*              nodalScalarValue;
  std::vector<fullVector<double> >* nodalVectorValue;

 public:
   Solution(const System& system);
  ~Solution(void);

  void write(const std::string name, Writer& writer) const;

 private:
  void interpolate(const FunctionSpace* fs);
};


/**
   @fn Solution::Solution
   @param system The System to use

   Instanciate a new Solution, based on the given System
   **

   @fn Solution::~Solution
   Deletes this Solution
   **
   
   @fn Solution::write
   @param name The file (@em without extension) where the 
   solution of the System will be written
   @param writer The Writer to use

   Writes the System's solution into the given file@n

   The file format is given by the Writer
   **
 */

#endif
