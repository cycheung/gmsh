#ifndef _PLOTBASIS_H_
#define _PLOTBASIS_H_

#include "GroupOfElement.h"
#include "fullMatrix.h"
#include "Basis.h"
#include "BasisScalar.h"
#include "BasisVector.h"

#include <fstream>
#include <string>
#include <vector>

class PlotBasis{
 private:
  const std::vector<MElement*>* element;
  const std::vector<MVertex*>*    node;
  
  int N;
  int E;
  int nFunction;

  std::vector<std::vector<double>*>*              nodalScalarValue;
  std::vector<std::vector<fullVector<double> >*>* nodalVectorValue;

  bool isScalar;

 public:
   PlotBasis(const GroupOfElement& group, const Basis& basis);
  ~PlotBasis(void);

  void write(const std::string name) const;

 private:
  void writeHeader(std::ofstream& out) const;
  void writeNodes(std::ofstream& out) const;
  void writeElements(std::ofstream& out) const;
  void writeNodalValues(std::ofstream& out, 
			const std::string name,
			int fun) const;

  void getGeometry(const GroupOfElement& group);
  void interpolate(const BasisScalar& basis);
  void interpolate(const BasisVector& basis);
};

#endif
