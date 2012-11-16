#include <cmath>

#include "Exception.h"
#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "Polynomial.h"
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
  gaussIntegration::get(goe.get(0).getType(), 2 * (order - 1), *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Function Space //
  fspace = new FunctionSpaceNode(goe, order);
}

FormulationLaplace::~FormulationLaplace(void){
  delete gC;
  delete gW;
  delete fspace;
}

double FormulationLaplace::weak(int dofI, int dofJ, 
				const GroupOfDof& god) const{
  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double> invJac(3, 3);        
  double integral = 0;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const vector<Polynomial>*> fun = 
    fspace->getGradLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = celement.getJacobian((*gC)(g, 0), 
				      (*gC)(g, 1), 
				      (*gC)(g, 2), 
				      invJac);
    invJac.invertInPlace();

    phiI = Mapper::grad(Polynomial::at(*fun[dofI], 
				       (*gC)(g, 0), 
				       (*gC)(g, 1),
				       (*gC)(g, 2)),
			invJac);

    phiJ = Mapper::grad(Polynomial::at(*fun[dofJ], 
				       (*gC)(g, 0), 
				       (*gC)(g, 1), 
			       (*gC)(g, 2)),
			invJac);
    
    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return integral;
}
