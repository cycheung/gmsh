#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include <string>
#include <vector>

#include "System.h"
#include "Writer.h"

#include "fullMatrix.h"
#include "MElement.h"
#include "DofManager.h"

class Solution{
 private:
  fullVector<double>* sol;

  const Mesh*                         mesh;
  const std::vector<const MElement*>* element;
  const DofManager*                   dofM;
  int                                 nDof;

  bool                              isScalar;
  std::vector<double>*              nodalScalarValue;
  std::vector<fullVector<double> >* nodalVectorValue;

 public:
   Solution(const System& system, const Mesh& mesh);
  ~Solution(void);

  void write(const std::string name, Writer& writer) const;

 private:
  void interpolateScalar(void);
  void interpolateVector(void);
};

#endif
