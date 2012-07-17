#include "Jacobian.h"
#include "Exception.h"

Jacobian::Jacobian(const MElement& element){
  this->element = const_cast<MElement*>(&element);
}

Jacobian::~Jacobian(void){
}

fullVector<double> Jacobian::map(const fullVector<double>& UVW) const{
  fullVector<double> XYZ(3);
  fullMatrix<double> jac(3, 3);

  element->getJacobian(UVW(0), UVW(1), UVW(2), jac);

  jac.mult(UVW, XYZ);

  //XY(0) = u * dxdu + v * dxdv + nodeX[0]; 
  //XY(1) = u * dydu + v * dydv + nodeY[0]; 

  return XYZ;
}

fullVector<double> Jacobian::grad(const fullVector<double>& gradUVW) const{
  fullVector<double> gradXYZ(3);
  fullMatrix<double> jac(3, 3);

  element->getJacobian(gradUVW(0), gradUVW(1), gradUVW(2), jac);
  
  jac.invertInPlace();
  jac.multWithATranspose(gradUVW, 1, 1, gradXYZ);
  
  //gradXY(0) = gradUV(0) * dudx + gradUV(1) * dvdx;
  //gradXY(1) = gradUV(0) * dudy + gradUV(1) * dvdy;    
  
  return gradXYZ;
}

fullVector<double> Jacobian::invMap(const fullVector<double>& XYZ) const{
  fullVector<double> UVW(3);
  
  throw Exception("invMap not implemented");

  //UV(0) = (XY(0) - nodeX[0]) * dudx + (XY(1) - nodeY[0]) * dudy;
  //UV(1) = (XY(0) - nodeX[0]) * dvdx + (XY(1) - nodeY[0]) * dvdy;  

  return UVW;
}
