#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "FormulationLaplace.h"
#include "Mapper.h"
#include "MElement.h"
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
				const GroupOfDof& god) const{

  fullMatrix<double>  invJac(3, 3);        
  MElement& element = const_cast<MElement&>(god.getGeoElement());
  double integral   = 0;

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = element.getJacobian((*gC)(g, 0), 
				     (*gC)(g, 1), 
				     (*gC)(g, 2), 
				     invJac);
    invJac.invertInPlace();

    fullVector<double> phiI = Mapper::grad(Polynomial::at(gradBasis[nodeI], 
							  (*gC)(g, 0), 
							  (*gC)(g, 1),
							  (*gC)(g, 2)),
					   invJac);
				       
    fullVector<double> phiJ = Mapper::grad(Polynomial::at(gradBasis[nodeJ], 
							  (*gC)(g, 0), 
							  (*gC)(g, 1), 
							  (*gC)(g, 2)),
					   invJac);

    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }
			
  return integral;
}
