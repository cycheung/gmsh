#include <cmath>

#include "GaussIntegration.h"
#include "Mapper.h"
#include "Polynomial.h"

#include "FormulationProjection.h"

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
  gaussIntegration::get(goe.get(0).getType(), 6, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Function Space //
  fspace = new FunctionSpaceEdge(goe, 3);
}

FormulationProjection::~FormulationProjection(void){
  delete gC;
  delete gW;
  delete fspace;
}

double FormulationProjection::weak(int edgeI, int edgeJ, 
				   const GroupOfDof& god) const{
  // Init //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double> invJac(3, 3);       
  double integral = 0;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const vector<Polynomial>*> fun = 
    fspace->getLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = celement.getJacobian((*gC)(g, 0), 
				      (*gC)(g, 1), 
				      (*gC)(g, 2), 
				      invJac);
    invJac.invertInPlace();

    phiI = Mapper::grad(Polynomial::at(*fun[edgeI],
				       (*gC)(g, 0), 
				       (*gC)(g, 1),
				       (*gC)(g, 2)),
			invJac);
    
    phiJ = Mapper::grad(Polynomial::at(*fun[edgeJ],
				       (*gC)(g, 0), 
				       (*gC)(g, 1),
				       (*gC)(g, 2)),
			invJac);

    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return integral;
}

double FormulationProjection::rhs(int equationI,
				  const GroupOfDof& god) const{
  // Init //
  fullMatrix<double>  invJac(3, 3);        
  fullVector<double>  phi(3);
  fullVector<double>& ff = const_cast<fullVector<double>&>(*f);
  double integral        = 0;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const vector<Polynomial>*> fun = 
    fspace->getLocalFunctions(element);  

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){  
    double det = celement.getJacobian((*gC)(g, 0), 
				      (*gC)(g, 1), 
				      (*gC)(g, 2), 
				      invJac);
    invJac.invertInPlace();

    phi = Mapper::grad(Polynomial::at(*fun[equationI],
				      (*gC)(g, 0), 
				      (*gC)(g, 1),
				      (*gC)(g, 2)),
		       invJac);
 
    integral += ff * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}
