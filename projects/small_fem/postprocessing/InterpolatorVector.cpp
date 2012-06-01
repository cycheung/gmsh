#include "InterpolatorVector.h"
#include "Exception.h"

InterpolatorVector::InterpolatorVector(void){
  scalar = false;
}

InterpolatorVector::~InterpolatorVector(void){
}

std::vector<fullVector<double>*>*
InterpolatorVector::getNodeValue(void) const{
  if(!gotInterpolation)
    throw Exception("Field not interpolated");

  return nodeValue;
}
