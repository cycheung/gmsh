#include "InterpolatorVector.h"
#include "Exception.h"

InterpolatorVector::InterpolatorVector(void){
  scalar           = false;
  gotInterpolation = false;
}

InterpolatorVector::~InterpolatorVector(void){
}

const std::vector<fullVector<double>*>&
InterpolatorVector::getNodeValue(void) const{
  if(!gotInterpolation)
    throw Exception("Field not interpolated");

  return *nodeValue;
}
