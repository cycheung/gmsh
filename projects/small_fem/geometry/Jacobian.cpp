#include "Jacobian.h"
#include "Exception.h"

Jacobian::Jacobian(const std::vector<Node*>& nodes){
  nNode = nodes.size();
  
  jac = new fullMatrix<double>(3, 3);

  nodeX = new double[nNode];
  nodeY = new double[nNode];
  nodeZ = new double[nNode];

  for(int i = 0; i < nNode; i++){
    nodeX[i] = nodes[i]->getX();
    nodeY[i] = nodes[i]->getY();
    nodeZ[i] = nodes[i]->getZ();
  }  

  switch(nNode){
  case 3: triJac(); break;

  default: 
    throw Exception
      ("I can't compute the Jacobian of an Element with %d nodes",
       nNode);
  }
}

Jacobian::~Jacobian(void){
  delete[] nodeX;
  delete[] nodeY;
  delete[] nodeZ;
  delete   jac;
}

fullVector<double> Jacobian::grad(const fullVector<double>& gradUV) const{
  fullVector<double> gradXY(2);
  
  gradXY(0) = gradUV(0) * dudx + gradUV(1) * dvdx;
  gradXY(1) = gradUV(0) * dudy + gradUV(1) * dvdy;    
  
  return gradXY;
}

fullVector<double> Jacobian::invMap(const fullVector<double>& XY) const{
  fullVector<double> UV(2);
  
  UV(0) = (XY(0) - nodeX[0]) * dudx + (XY(1) - nodeY[0]) * dudy;
  UV(1) = (XY(0) - nodeX[0]) * dvdx + (XY(1) - nodeY[0]) * dvdy;  
  
  return UV;
}

fullVector<double> Jacobian::invMap(const double x, const double y) const{
  fullVector<double> UV(2);
  
  UV(0) = (x - nodeX[0]) * dudx + (y - nodeY[0]) * dudy;
  UV(1) = (x - nodeX[0]) * dvdx + (y - nodeY[0]) * dvdy;  
  
  return UV;
}

fullVector<double> Jacobian::map(const fullVector<double>& UV) const{
  fullVector<double> XY(2);

  XY(0) = UV(0) * dxdu + UV(1) * dxdv + nodeX[0];
  XY(1) = UV(0) * dydu + UV(1) * dydv + nodeY[0];  

  return XY;
}

fullVector<double> Jacobian::map(const double u, const double v) const{
  fullVector<double> XY(2);

  XY(0) = u * dxdu + v * dxdv + nodeX[0];
  XY(1) = u * dydu + v * dydv + nodeY[0];  

  return XY;
}

void Jacobian::triJac(void){

  dxdu = nodeX[1] - nodeX[0];
  dxdv = nodeX[2] - nodeX[0];
  dydu = nodeY[1] - nodeY[0];
  dydv = nodeY[2] - nodeY[0];
  
  detDxDu = (dxdu * dydv) - (dxdv * dydu);
  
  dudx = +dydv / detDxDu;
  dudy = -dxdv / detDxDu;
  dvdx = -dydu / detDxDu;
  dvdy = +dxdu / detDxDu; 
}
