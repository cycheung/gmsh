#include "Mapper.h"
#include "Exception.h"

Mapper::Mapper(void){
}

Mapper::~Mapper(void){
}

fullVector<double> Mapper::invMap(const fullVector<double>& XYZ, 
				  const fullMatrix<double>& invJac){
  fullVector<double> UVW(3);
  
  throw Exception("invMap not implemented");

  //UV(0) = (XY(0) - nodeX[0]) * dudx + (XY(1) - nodeY[0]) * dudy;
  //UV(1) = (XY(0) - nodeX[0]) * dvdx + (XY(1) - nodeY[0]) * dvdy;  

  return UVW;
}
