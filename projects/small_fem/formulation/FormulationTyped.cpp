#include <complex>
#include "FormulationTyped.h"

template<>
std::string FormulationTyped<double>::getType(void) const{
  return "real";
}

template<>
std::string FormulationTyped<std::complex<double> >::getType(void) const{
  return "complex";
}
