#include "SystemAbstract.h"

template<>
const double SystemAbstract<double>::minusSign = -1;

template<>
const std::complex<double> SystemAbstract<std::complex<double> >::minusSign =
  std::complex<double>(-1, 0);
