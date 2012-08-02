#ifndef _WRITER_H_
#define _WRITER_H_

#include <string>

class Writer{
 public:
  Writer(void);

  virtual ~Writer(void);

  virtual void write(const std::string name) const = 0;
};

#endif
