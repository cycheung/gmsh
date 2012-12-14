#include "Mapper.h"

Mapper::Mapper(void){
}

Mapper::~Mapper(void){
}

// WARNING                //
// invJac is Transposed ! //

fullVector<double> Mapper::map(const fullVector<double>& UVW,
			       const fullVector<double>& originXYZ,
			       const fullMatrix<double>& jac){
  fullVector<double> XYZ(3);

  jac.multWithATranspose(UVW, 1, 0, XYZ);
  XYZ.axpy(originXYZ, +1);

  return XYZ;
}

fullVector<double> Mapper::invMap(const fullVector<double>& XYZ, 
				  const fullVector<double>& originXYZ,
				  const fullMatrix<double>& invJac){
  fullVector<double> UVW(3);
  fullVector<double> sub(XYZ);
  sub.axpy(originXYZ, -1);

  invJac.multWithATranspose(sub, 1, 0, UVW);
  
  return UVW;
}

fullVector<double> Mapper::grad(const fullVector<double>& gradUVW, 
				const fullMatrix<double>& invJac){
  fullVector<double> gradXYZ(3);
  
  invJac.mult(gradUVW, gradXYZ);
  return gradXYZ;
}

fullVector<double> Mapper::grad(double gradU,
				double gradV,
				double gradW, 
				const fullMatrix<double>& invJac){
  
  fullVector<double> gradXYZ(3);
  fullVector<double> gradUVW(3);
  
  gradUVW(0) = gradU;
  gradUVW(1) = gradV;
  gradUVW(2) = gradW;

  invJac.mult(gradUVW, gradXYZ);
  return gradXYZ;
}

fullVector<double> Mapper::curl(const fullVector<double>& curlUVW, 
				const fullMatrix<double>& jac,
				double invDet){
  fullVector<double> curlXYZ(3);
  
  jac.multWithATranspose(curlUVW, invDet, 0, curlXYZ);
  
  return curlXYZ;
}

fullVector<double> Mapper::curl(double curlU,
				double curlV,
				double curlW, 
				const fullMatrix<double>& jac,
				double invDet){
  fullVector<double> curlXYZ(3);
  fullVector<double> curlUVW(3);
  
  curlUVW(0) = curlU;
  curlUVW(1) = curlV;
  curlUVW(2) = curlW;
  
  jac.multWithATranspose(curlUVW, invDet, 0, curlXYZ);
  
  return curlXYZ;
}
