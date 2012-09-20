#include <cmath>
#include "PoissonSquare.h"

unsigned int PoissonSquare::max = 225;

double PoissonSquare::pi =
  3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982;

PoissonSquare::PoissonSquare(const GroupOfElement& goe){
  scalar = true;
  domain = &goe; 
  
  compute();
}

PoissonSquare::~PoissonSquare(void){
}

double PoissonSquare::fScalar(double x, double y, double z){
  double res = (1 - (x * x)) / 2;  

  for(unsigned int k = 1; k <= max; k += 2){
    res -= 
      16 / (pi * pi * pi) * 
      
      ((sin(k * pi  * (1 + x) / 2)) / 
       (k * k * k * sinh(k * pi))) *
      
      (sinh(k * pi * (1 + y) / 2) + sinh(k * pi * (1 - y) / 2));
  }
  

  return res;
}
