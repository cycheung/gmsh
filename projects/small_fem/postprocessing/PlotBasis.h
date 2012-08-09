#ifndef _PLOTBASIS_H_
#define _PLOTBASIS_H_

#include "Writer.h"

#include "GroupOfElement.h"
#include "fullMatrix.h"
#include "Basis.h"

#include <string>
#include <vector>

class PlotBasis{
 private:
  Writer* writer;
  int nFunction;
  bool isScalar;

  int N;
  int E;

  const std::vector<const MElement*>* element;
  const std::vector<MVertex*>*        node;

  std::vector<double>**              nodalScalarValue;
  std::vector<fullVector<double> >** nodalVectorValue;

 public:
  PlotBasis(const GroupOfElement& group, 
	    const Basis& basis,
	    Writer& writer);
  
  virtual ~PlotBasis(void);

  virtual void plot(const std::string name) const;

 private:
  void getGeometry(const GroupOfElement& group);
  void interpolate(const BasisScalar& basis);
  void interpolate(const BasisVector& basis);
};

#endif
