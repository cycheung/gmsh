#ifndef _WRITER_H_
#define _WRITER_H_

#include "fullMatrix.h"
#include "MElement.h"
#include "MVertex.h"

#include <string>
#include <vector>

class Writer{
 protected:
  bool hasDomain;
  bool hasValue;
  bool isScalar;

  int  N;
  int  E;

  const std::vector<const MElement*>* element;
  const std::vector<MVertex*>*        node;

  std::vector<double>*              nodalScalarValue;
  std::vector<fullVector<double> >* nodalVectorValue;

 public:
  Writer(void);

  virtual ~Writer(void);

  virtual void write(const std::string name) const = 0;

  void setValues(std::vector<double>& value);
  void setValues(std::vector<fullVector<double> >& value);
  void setDomain(const std::vector<const MElement*>& element);
};

#endif
