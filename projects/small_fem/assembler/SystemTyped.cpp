#include <complex>
#include "SystemTyped.h"

template<>
std::string SystemTyped<double>::
getType(void) const{
  return "real";
}

template<>
std::string SystemTyped<std::complex<double> >::
getType(void) const{
  return "complex";
}
