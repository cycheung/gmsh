#include <cmath>

#include "Exception.h"
#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "BasisScalar.h"
#include "Mapper.h"

#include "FormulationLaplace.h"

using namespace std;

FormulationLaplace::FormulationLaplace(const GroupOfElement& goe,
				       unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw 
      Exception("Can't have a Laplace formulation of order 0");

  // Gaussian Quadrature Data //
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), order, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Function Space //
  fspace = new FunctionSpaceNode(goe, order);
}

FormulationLaplace::~FormulationLaplace(void){
  delete gC;
  delete gW;
  delete fspace;
}

double FormulationLaplace::weak(int nodeI, int nodeJ, 
				const GroupOfDof& god) const{
  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double> invJac(3, 3);        
  double integral = 0;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const Polynomial*> fun = fspace->getLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = celement.getJacobian((*gC)(g, 0), 
				      (*gC)(g, 1), 
				      (*gC)(g, 2), 
				      invJac);
    invJac.invertInPlace();

    phiI = Mapper::grad(Polynomial::at(fun[nodeI]->gradient(), 
				       (*gC)(g, 0), 
				       (*gC)(g, 1),
				       (*gC)(g, 2)),
			invJac);

    phiJ = Mapper::grad(Polynomial::at(fun[nodeJ]->gradient(), 
				       (*gC)(g, 0), 
				       (*gC)(g, 1), 
			       (*gC)(g, 2)),
			invJac);
    
    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return integral;
  /*
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
  */
}
