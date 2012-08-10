#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include <string>
#include <vector>

#include "System.h"
#include "Writer.h"

#include "fullMatrix.h"
#include "MElement.h"
#include "DofManager.h"

#include "FunctionSpace.h"
#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

class Solution{
 private:
  fullVector<double>* sol;

  const Mesh*                         mesh;
  const std::vector<const MElement*>* element;
  const std::vector<GroupOfDof*>*     god;
  const DofManager*                   dofM;
  const FunctionSpace*                fs;

  bool                              isScalar;
  std::vector<double>*              nodalScalarValue;
  std::vector<fullVector<double> >* nodalVectorValue;

 public:
   Solution(const System& system, const Mesh& mesh);
  ~Solution(void);

  void write(const std::string name, Writer& writer) const;

 private:
  void interpolateScalar(const FunctionSpaceScalar& fs);
  void interpolateVector(const FunctionSpaceVector& fs);
};

#endif
