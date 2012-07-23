#include "InterpolatorScalar.h"
#include "Exception.h"

InterpolatorScalar::InterpolatorScalar(void){
  scalar           = true;
  gotInterpolation = false;
}

InterpolatorScalar::~InterpolatorScalar(void){
}

const std::vector<double>&
InterpolatorScalar::getNodeValue(void) const{
  if(!gotInterpolation)
    throw Exception("Field not interpolated");

  return *nodeValue;
}
