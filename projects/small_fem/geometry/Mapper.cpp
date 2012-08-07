#include "Mapper.h"
#include "Exception.h"

Mapper::Mapper(void){
}

Mapper::~Mapper(void){
}

fullVector<double> Mapper::map(const fullVector<double>& UVW,
			       const fullMatrix<double>& jac){
  // WARNING             //
  // jac is Transposed ! //

  fullVector<double> XYZ(3);

  throw Exception("Bad implementation of Mapper::map");

  //jac.multWithATranspose(UVW, 1, 0, XYZ);
  // + origin !!!!!!!!
  return XYZ;
}

fullVector<double> Mapper::invMap(const fullVector<double>& XYZ, 
				  const fullVector<double>& origin,
				  const fullMatrix<double>& invJac){
  // WARNING                //
  // invJac is Transposed ! //

  fullVector<double> UVW(3);
  fullVector<double> sub(XYZ);
  sub.axpy(origin, -1);

  invJac.multWithATranspose(sub, 1, 0, UVW);
  
  return UVW;
}

fullVector<double> Mapper::grad(const fullVector<double>& gradUVW, 
				const fullMatrix<double>& invJac){
  // WARNING                //
  // invJac is Transposed ! //

  fullVector<double> gradXYZ(3);
  
  invJac.mult(gradUVW, gradXYZ);
  return gradXYZ;
}
