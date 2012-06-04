#include "fullMatrix.h"
#include "FormulationLaplace.h"
#include <cmath>

using namespace std;

FormulationLaplace::FormulationLaplace(void){
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
  // Generate Basis
  base = new TriNodeBasis(1);  
  const Polynomial* basis = base->getBasis();

  // Take gradient
  basisSize = base->getSize();  
  gradBasis = new vector<Polynomial> [basisSize];

  for(int i = 0; i < basisSize; i++)
    gradBasis[i] = basis[i].gradient();

  // Interpolator //
  interp = new InterpolatorNode();
}

FormulationLaplace::~FormulationLaplace(void){
  delete   interp;
  delete   base;
  delete[] gradBasis;
}

double FormulationLaplace::weak(const int nodeI, const int nodeJ, 
				const GroupOfDof& god) const{
  const Jacobian& jac = god.getJacobian();

  double integral = 0;  
  for(int g = 0; g < G; g++){
    fullVector<double> phiI = jac.grad(Polynomial::at(gradBasis[nodeI], gx[g], gy[g], 0));
    fullVector<double> phiJ = jac.grad(Polynomial::at(gradBasis[nodeJ], gx[g], gy[g], 0));

    integral += phiI * phiJ * fabs(jac.det()) * gw[g];
  }

  return integral;
}
