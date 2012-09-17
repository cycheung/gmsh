#include <cmath>

#include "Exception.h"
#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "Polynomial.h"
#include "Mapper.h"

#include "FormulationPoisson.h"

using namespace std;

FormulationPoisson::FormulationPoisson(const GroupOfElement& goe,
				       unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw 
      Exception("Can't have a Poisson formulation of order 0");

  // Gaussian Quadrature Data //
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), order + 4, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points
  /*
  cout << "Gauss: " << G << endl;
  
  for(int i = 0; i < G; i++)
    cout << "  -- [" 
	 << (*gC)(i, 0) << ", "
	 << (*gC)(i, 1) << ", "
	 << (*gC)(i, 2) << "]" << endl;
  */
  // Function Space //
  fspace = new FunctionSpaceNode(goe, order);
}

FormulationPoisson::~FormulationPoisson(void){
  delete gC;
  delete gW;
  delete fspace;
}

double FormulationPoisson::weak(int nodeI, int nodeJ, 
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
}

double FormulationPoisson::rhs(int equationI,
			       const GroupOfDof& god) const{

  // Init Some Stuff //
  fullMatrix<double> jac(3, 3);        
  double integral = 0;
  double phi;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const Polynomial*> fun = fspace->getLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = celement.getJacobian((*gC)(g, 0), 
				      (*gC)(g, 1), 
				      (*gC)(g, 2), 
				      jac);

    phi = fun[equationI]->at((*gC)(g, 0), 
			     (*gC)(g, 1), 
			     (*gC)(g, 2));

    integral += 1 * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}
