#include "Jacobian.h"

Jacobian::Jacobian(const MElement& element){
  this->element = const_cast<MElement*>(&element);

  //!@todo 
  //!const_cast is dirty !!
  //!Change MElement with const qualifier
}

Jacobian::~Jacobian(void){
}

double Jacobian::det(const fullVector<double>& UV) const{
  return 42;
}

double Jacobian::det(void) const{
  return 42;
}

fullVector<double> Jacobian::map(const fullVector<double>& UV) const{
  fullVector<double> XY(3);
  fullMatrix<double> jac(3, 3);

  element->getJacobian(UV(0), UV(1), UV(2), jac);

  jac.mult(UV, XY);

  return XY;
}

fullVector<double> Jacobian::grad(const fullVector<double>& gradUV) const{
  fullVector<double> gradXY(3);
  /*  
  gradXY(0) = gradUV(0) * dudx + gradUV(1) * dvdx;
  gradXY(1) = gradUV(0) * dudy + gradUV(1) * dvdy;    
  */
  return gradXY;
}

fullVector<double> Jacobian::invMap(const fullVector<double>& XY) const{
  fullVector<double> UV(3);
  /*
  UV(0) = (XY(0) - nodeX[0]) * dudx + (XY(1) - nodeY[0]) * dudy;
  UV(1) = (XY(0) - nodeX[0]) * dvdx + (XY(1) - nodeY[0]) * dvdy;  
  */
  return UV;
}
