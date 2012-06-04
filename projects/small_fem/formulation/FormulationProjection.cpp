#include "FormulationProjection.h"
#include <cmath>

using namespace std;

FormulationProjection::FormulationProjection(fullVector<double>& vectorToProject){
  // Vector to Project //
  f = &vectorToProject;

  // Gaussian Quadrature Data //
  G     = 4;

  gx[0] = 0.333333333333333;
  gx[1] = 0.6;
  gx[2] = 0.2;
  gx[3] = 0.2;

  gy[0] = 0.333333333333333;
  gy[1] = 0.2;
  gy[2] = 0.6;
  gy[3] = 0.2;

  gw[0] = -0.28125;
  gw[1] = +0.260416666666;
  gw[2] = +0.260416666666;
  gw[3] = +0.260416666666;

  // Basis //
  baseGen   = new TriNedelecBasis;
  basis     = baseGen->getBasis();
  basisSize = baseGen->getSize(); 

  // Interpolator //
  interp = new InterpolatorEdge(*baseGen);
}

FormulationProjection::~FormulationProjection(void){
  delete baseGen;
  delete interp;
}

double FormulationProjection::weak(const int edgeI, const int edgeJ, 
				   const GroupOfDof& god) const{
 
  const Jacobian& jac = god.getJacobian();
  int orientationI    = god.getOrientation(edgeI);
  int orientationJ    = god.getOrientation(edgeJ);
  int orientation     = orientationI * orientationJ;
  
  // Loop over Integration Point //
  double integral = 0;  
  for(int g = 0; g < G; g++){
    fullVector<double> phiI = jac.grad(Polynomial::at(basis[edgeI], gx[g], gy[g], 0));
    fullVector<double> phiJ = jac.grad(Polynomial::at(basis[edgeJ], gx[g], gy[g], 0));

    integral += phiI * phiJ * fabs(jac.det()) * gw[g] * orientation;
  }

  return integral;
}

double FormulationProjection::rhs(const int equationI,
				  const GroupOfDof& god) const{
 
  const Jacobian& jac = god.getJacobian();
  int orientation = god.getOrientation(equationI);

  // Loop over Integration Point //
  double integral = 0;
  for(int g = 0; g < G; g++){  
    fullVector<double> jPhiI = jac.grad(Polynomial::at(basis[equationI], gx[g], gy[g], 0));
 
    integral += (*f) * jPhiI * fabs(jac.det()) * gw[g] * orientation;
  }

  return integral;
}
