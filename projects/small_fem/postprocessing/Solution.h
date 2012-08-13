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

#endif
