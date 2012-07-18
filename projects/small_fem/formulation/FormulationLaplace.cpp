#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "FormulationLaplace.h"
#include <cmath>

using namespace std;

FormulationLaplace::FormulationLaplace(void){
  // Gaussian Quadrature Data //
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  gaussIntegration::getTriangle(1, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Basis //
  // Generate Basis
  base = new TriNodeBasis(1);  
  const vector<Polynomial>& basis = base->getBasis();

  // Take gradient
  basisSize = base->getSize();  
  gradBasis = new vector<Polynomial> [basisSize];

  for(int i = 0; i < basisSize; i++)
    gradBasis[i] = basis[i].gradient();

  // Interpolator //
  //interp = new InterpolatorNode();
}

FormulationLaplace::~FormulationLaplace(void){
  delete   gC;
  delete   gW;
  //delete   interp;
  delete   base;
  delete[] gradBasis;
}

double FormulationLaplace::weak(const int nodeI, const int nodeJ, 
				const GeoDof& god) const{

  double integral = 0;
  for(int g = 0; g < G; g++){
    fullVector<double> phiI = god.grad(Polynomial::at(gradBasis[nodeI], 
						      (*gC)(g, 0), 
						      (*gC)(g, 1),
						      (*gC)(g, 2)));
				       
    fullVector<double> phiJ = god.grad(Polynomial::at(gradBasis[nodeJ], 
						      (*gC)(g, 0), 
						      (*gC)(g, 1), 
						      (*gC)(g, 2)));

    integral += 
      phiI * phiJ * 
      fabs(god.getJacobian((*gC)(g, 0), 
			   (*gC)(g, 1), 
			   (*gC)(g, 2))) * (*gW)(g);
  }
			
  return integral;
}
