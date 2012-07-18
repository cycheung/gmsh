#include "GaussIntegration.h"
#include "FormulationProjection.h"
#include "Mapper.h"
#include "MElement.h"
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
				   const GroupOfDof& god) const{

  fullMatrix<double>  invJac(3, 3);        
  MElement& element = const_cast<MElement&>(god.getElement());
  double integral   = 0;
  int orientation   = 
    god.getOrientation(edgeI) * 
    god.getOrientation(edgeJ);
  
  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = element.getJacobian((*gC)(g, 0), 
				     (*gC)(g, 1), 
				     (*gC)(g, 2), 
				     invJac);
    invJac.invertInPlace();

    fullVector<double> phiI = Mapper::grad(Polynomial::at((*basis)[edgeI],
							  (*gC)(g, 0), 
							  (*gC)(g, 1),
							  (*gC)(g, 2)),
					   invJac);
    
    fullVector<double> phiJ = Mapper::grad(Polynomial::at((*basis)[edgeJ],
							  (*gC)(g, 0), 
							  (*gC)(g, 1),
							  (*gC)(g, 2)),
					   invJac);

    integral += phiI * phiJ * fabs(det) * (*gW)(g) * orientation;
  }

  return integral;
}

double FormulationProjection::rhs(const int equationI,
				  const GroupOfDof& god) const{
 
  fullMatrix<double>  invJac(3, 3);        
  MElement& element = const_cast<MElement&>(god.getElement());
  int orientation   = god.getOrientation(equationI);
  double integral   = 0;

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){  
    double det = element.getJacobian((*gC)(g, 0), 
				     (*gC)(g, 1), 
				     (*gC)(g, 2), 
				     invJac);
    invJac.invertInPlace();

    fullVector<double> phiI = Mapper::grad(Polynomial::at((*basis)[equationI],
							  (*gC)(g, 0), 
							  (*gC)(g, 1),
							  (*gC)(g, 2)),
					   invJac);
 
    integral += (*f) * phiI * fabs(det) * (*gW)(g) * orientation;
  }

  return integral;
}
