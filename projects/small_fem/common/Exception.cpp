#include <cstdio>
#include "Exception.h"

Exception::Exception(const std::string format, ...){
  char    buf[256];
  va_list list;
  
  va_start(list, format);
  
  if(vsnprintf(buf, 256, format.c_str(), list) < 0)
    why = new std::string("Unknown Exception");

  else
    why = new std::string(buf);
}

Exception::~Exception(void) throw(){
  delete why;
}
  
const char* Exception::what(void) const throw(){
  return why->c_str();
}
