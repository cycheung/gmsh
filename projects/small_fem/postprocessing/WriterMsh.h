#ifndef _WRITERMSH_H_
#define _WRITERMSH_H_

#include "Writer.h"

#include "MElement.h"
#include "MVertex.h"
#include "fullMatrix.h"

#include <vector>
#include <fstream>
#include <string>

class WriterMsh: public Writer{
 protected:
  int  N;
  int  E;
  bool isScalar;

  const std::vector<MElement*>* element;
  const std::vector<MVertex*>*    node;
  
  std::vector<double>**              nodalScalarValue;
  std::vector<fullVector<double> >** nodalVectorValue;
    
  mutable std::ofstream* out;
  
 public:
  WriterMsh(void);
  
  virtual ~WriterMsh(void);
  
  virtual void write(const std::string name) const = 0;

 protected:
  void writeHeader(void) const;
  void writeNodes(void) const;
  void writeElements(void) const;
  void writeNodalValues(const std::string name, int num) const;
};

#endif
