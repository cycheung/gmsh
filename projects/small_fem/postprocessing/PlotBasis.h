#ifndef _PLOTBASIS_H_
#define _PLOTBASIS_H_

#include "WriterMsh.h"

#include "GroupOfElement.h"
#include "fullMatrix.h"
#include "Basis.h"
#include "BasisScalar.h"
#include "BasisVector.h"

#include <fstream>
#include <string>
#include <vector>

class PlotBasis: public WriterMsh{
 private:
  int nFunction;

 public:
  PlotBasis(const GroupOfElement& group, const Basis& basis);
  
  virtual ~PlotBasis(void);

  virtual void write(const std::string name) const;

 private:
  void getGeometry(const GroupOfElement& group);
  void interpolate(const BasisScalar& basis);
  void interpolate(const BasisVector& basis);
};

#endif
