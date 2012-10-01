#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include <string>
#include <vector>

#include "System.h"
#include "Writer.h"
#include "fullMatrix.h"

#include "Mesh.h"
#include "GroupOfElement.h"
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

  double (*fScalar)(fullVector<double>& xyz);

  int fsType;
  
  bool                              scalar;
  std::vector<double>*              nodalScalarValue;
  std::vector<fullVector<double> >* nodalVectorValue;
  const GroupOfElement*             visuDomain;

 public:
   Solution(const System& system);
   Solution(const System& system, const GroupOfElement& visu);
   Solution(double (*f)(fullVector<double>& xyz), const GroupOfElement& visu);
  ~Solution(void);

  void write(const std::string name, Writer& writer) const;
  bool isScalar(void) const;
  
  std::vector<double>&              getNodalScalarValue(void) const;
  std::vector<fullVector<double> >& getNodalVectorValue(void) const;

 private:
  void init(const System& system);
  void interpolate(void);
  void interpolateOnVisu(void); 
  void evaluateF(void);
};


/**
   @fn Solution::Solution
   @param system The System to use

   Instanciate a new Solution, based on the given System

   @note
   The interpolation will be done on the  @em support of the
   FunctionSpace.
   **

   @fn Solution::Solution
   @param system The System to use
   @param visu The GroupOfElement to use for interpolation

   Instanciate a new Solution, based on the given System

   @note
   The interpolation will be done on the given GroupOfElement
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

//////////////////////
// Inline Functions //
//////////////////////

inline bool Solution::isScalar(void) const{
  return scalar;
}

#endif
