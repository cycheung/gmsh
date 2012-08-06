#ifndef _WRITERMSH_H_
#define _WRITERMSH_H_

#include "Writer.h"
#include <fstream>

class WriterMsh: public Writer{
 protected:
  mutable std::ofstream* out;
  
 public:
  WriterMsh(void);
  
  virtual ~WriterMsh(void);
  
  virtual void write(const std::string name) const;

 protected:
  void writeHeader(void) const;
  void writeNodes(void) const;
  void writeElements(void) const;
  void writeNodalValues(const std::string name) const;
};

#endif
