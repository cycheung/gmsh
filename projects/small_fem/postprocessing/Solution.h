#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include <string>
#include <vector>

#include "System.h"
#include "Writer.h"

#include "fullMatrix.h"
#include "DofManager.h"
#include "MElement.h"
#include "GroupOfDof.h"

#include "FunctionSpace.h"
#include "FunctionSpaceScalar.h"

class Solution{
 private:
  fullVector<double>* sol;

  const DofManager*               dofM;
  const std::vector<GroupOfDof*>* god;
  const std::vector<MElement*>*   element;
  int                             nGod;
  int                             nVertex;

  bool                              isScalar;
  const FunctionSpaceScalar*        fsScalar;

  std::vector<double>*              nodalScalarValue;
  std::vector<fullVector<double> >* nodalVectorValue;

 public:
   Solution(const System& system);
  ~Solution(void);

  void write(const std::string name, Writer& writer) const;

 private:
  void interpolateScalar(void);
  void interpolateVector(void);
};

#endif
