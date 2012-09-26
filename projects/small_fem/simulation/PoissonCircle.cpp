#include "PoissonCircle.h"

PoissonCircle::PoissonCircle(const GroupOfElement& goe){
  scalar = true;
  domain = &goe; 
  
  compute();
}

PoissonCircle::~PoissonCircle(void){
}

double PoissonCircle::fScalar(double x, double y, double z){
  return (x * x + y * y - 1) / 4;
}
