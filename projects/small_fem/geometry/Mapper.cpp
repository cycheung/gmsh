#include "Mapper.h"
#include "Exception.h"

Mapper::Mapper(void){
}

Mapper::~Mapper(void){
}

fullVector<double> Mapper::map(const fullVector<double>& UVW,
				 MElement& element){
  fullVector<double> XYZ(3);
  fullMatrix<double> jac(3, 3);

  element.getJacobian(UVW(0), UVW(1), UVW(2), jac);

  jac.mult(UVW, XYZ);

  //XY(0) = u * dxdu + v * dxdv + nodeX[0]; 
  //XY(1) = u * dydu + v * dydv + nodeY[0]; 

  return XYZ;
}

fullVector<double> Mapper::grad(const fullVector<double>& gradUVW, 
				  MElement& element){
  fullVector<double> gradXYZ(3);
  fullMatrix<double> jac(3, 3);

  element.getJacobian(gradUVW(0), gradUVW(1), gradUVW(2), jac);
  
  jac.invertInPlace();
  jac.multWithATranspose(gradUVW, 1, 1, gradXYZ);
  
  //gradXY(0) = gradUV(0) * dudx + gradUV(1) * dvdx;
  //gradXY(1) = gradUV(0) * dudy + gradUV(1) * dvdy;    
  
  return gradXYZ;
}

fullVector<double> Mapper::invMap(const fullVector<double>& XYZ, 
				    MElement& element){
  fullVector<double> UVW(3);
  
  throw Exception("invMap not implemented");

  //UV(0) = (XY(0) - nodeX[0]) * dudx + (XY(1) - nodeY[0]) * dudy;
  //UV(1) = (XY(0) - nodeX[0]) * dvdx + (XY(1) - nodeY[0]) * dvdy;  

  return UVW;
}
