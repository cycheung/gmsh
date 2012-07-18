#include "GaussIntegration.h"
#include "FormulationProjection.h"
#include <cmath>

using namespace std;

FormulationProjection::FormulationProjection(fullVector<double>& vectorToProject){
  // Vector to Project //
  f = &vectorToProject;

  // Gaussian Quadrature Data //
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  gaussIntegration::getTriangle(2, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Basis //
  baseGen   = new TriNedelecBasis;
  basis     = &(baseGen->getBasis());
  basisSize = baseGen->getSize(); 

  // Interpolator //
  //interp = new InterpolatorEdge(*baseGen);
}

FormulationProjection::~FormulationProjection(void){
  delete gC;
  delete gW;
  delete baseGen;
  //delete interp;
}

double FormulationProjection::weak(const int edgeI, const int edgeJ, 
				   const GeoDof& god) const{

  int orientation = 
    god.getOrientation(edgeI) * 
    god.getOrientation(edgeJ);
  
  // Loop over Integration Point //
  double integral = 0;  
  for(int g = 0; g < G; g++){
    fullVector<double> phiI = god.grad(Polynomial::at((*basis)[edgeI],
						      (*gC)(g, 0), 
						      (*gC)(g, 1),
						      (*gC)(g, 2)));
    
    fullVector<double> phiJ = god.grad(Polynomial::at((*basis)[edgeJ],
						      (*gC)(g, 0), 
						      (*gC)(g, 1),
						      (*gC)(g, 2)));

    integral += 
      phiI * phiJ * 
      fabs(god.getJacobian((*gC)(g, 0), 
			   (*gC)(g, 1), 
			   (*gC)(g, 2))) * (*gW)(g) * orientation;
  }

  return integral;
}

double FormulationProjection::rhs(const int equationI,
				  const GeoDof& god) const{
 
  int orientation = god.getOrientation(equationI);

  // Loop over Integration Point //
  double integral = 0;
  for(int g = 0; g < G; g++){  
    fullVector<double> phiI = god.grad(Polynomial::at((*basis)[equationI],
						      (*gC)(g, 0), 
						      (*gC)(g, 1),
						      (*gC)(g, 2)));
 
    integral += 
      (*f) * phiI *
      fabs(god.getJacobian((*gC)(g, 0), 
			   (*gC)(g, 1), 
			   (*gC)(g, 2))) * (*gW)(g) * orientation;
  }

  return integral;
}
