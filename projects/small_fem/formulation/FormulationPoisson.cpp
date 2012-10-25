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

  // Gaussian Quadrature Data (LHS) // 
  // NB: We need to integrad a grad * grad !
  //     and order(grad f) = order(f) - 1
  gCL = new fullMatrix<double>();
  gWL = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), (order - 1) + (order - 1) , *gCL, *gWL);

  GL = gWL->size(); // Nbr of Gauss points

  // Gaussian Quadrature Data (RHS) //
  // NB: We need to integrad a f * g !
  //     and here, g = 1
  gCR = new fullMatrix<double>();
  gWR = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), order, *gCR, *gWR);

  GR = gWR->size(); // Nbr of Gauss points


  // Function Space //
  fspace = new FunctionSpaceNode(goe, order);
}

FormulationPoisson::~FormulationPoisson(void){
  delete gCL;
  delete gWL;
  delete gCR;
  delete gWR;
  delete fspace;
}

double FormulationPoisson::weak(int dofI, int dofJ, 
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
  for(int g = 0; g < GL; g++){
    double det = celement.getJacobian((*gCL)(g, 0), 
				      (*gCL)(g, 1), 
				      (*gCL)(g, 2), 
				      invJac);
    invJac.invertInPlace();

    phiI = Mapper::grad(Polynomial::at(*fun[dofI], 
				       (*gCL)(g, 0), 
				       (*gCL)(g, 1),
				       (*gCL)(g, 2)),
			invJac);

    phiJ = Mapper::grad(Polynomial::at(*fun[dofJ], 
				       (*gCL)(g, 0), 
				       (*gCL)(g, 1), 
				       (*gCL)(g, 2)),
			invJac);
    
    integral += phiI * phiJ * fabs(det) * (*gWL)(g);
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
  for(int g = 0; g < GR; g++){
    double det = celement.getJacobian((*gCR)(g, 0), 
				      (*gCR)(g, 1), 
				      (*gCR)(g, 2), 
				      jac);

    phi = fun[equationI]->at((*gCR)(g, 0), 
			     (*gCR)(g, 1), 
			     (*gCR)(g, 2));

    integral += -1 * phi * fabs(det) * (*gWR)(g);
  }

  return integral;
}
