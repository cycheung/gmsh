#include "GaussIntegration.h"
#include "Mapper.h"
#include "BasisVector.h"

#include "FunctionSpaceEdge.h"
#include "FormulationProjection.h"

#include <cmath>

using namespace std;

FormulationProjection::FormulationProjection(const GroupOfElement& goe,
					     const fullVector<double>& vectorToProject){
  // Vector to Project //
  f = &vectorToProject;

  // Gaussian Quadrature Data //
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), 2, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Function Space //
  fspace = new FunctionSpaceEdge(goe, 0);

  // Basis //
  const BasisVector& base = 
    static_cast<const BasisVector&>(fspace->getBasis(goe.get(0)));
  basis = &(base.getBasis());

  basisSize = base.getSize(); 
}

FormulationProjection::~FormulationProjection(void){
  delete gC;
  delete gW;
  delete fspace;
}

double FormulationProjection::weak(const int edgeI, const int edgeJ, 
				   const GroupOfDof& god) const{

  fullMatrix<double>  invJac(3, 3);        
  MElement& element = const_cast<MElement&>(god.getGeoElement());

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
  MElement& element      = const_cast<MElement&>(god.getGeoElement());
  fullVector<double>& ff = const_cast<fullVector<double>&>(*f);

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
 
    integral += ff * phiI * fabs(det) * (*gW)(g) * orientation;
  }

  return integral;
}
