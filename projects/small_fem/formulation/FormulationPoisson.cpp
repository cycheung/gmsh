#include <cmath>

#include "Exception.h"
#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "BasisScalar.h"
#include "Mapper.h"

#include "FunctionSpaceNode.h"
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
  gaussIntegration::get(goe.get(0).getType(), order, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Function Space //
  FunctionSpaceNode* fspace = new FunctionSpaceNode(goe, order);
  this->fspace              = fspace;

  // Basis //
  // Get Basis
  const BasisScalar& base = fspace->getBasis(goe.get(0));
  basis = &(base.getFunctions());

  // Take gradient
  unsigned int basisSize = basis->size();

  gradBasis = new vector<Polynomial>[basisSize];

  for(unsigned int i = 0; i < basisSize; i++)
    gradBasis[i] = (*basis)[i].gradient();
}

FormulationPoisson::~FormulationPoisson(void){
  delete   gC;
  delete   gW;
  delete   fspace;
  delete[] gradBasis;
}

double FormulationPoisson::weak(int nodeI, int nodeJ, 
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

double FormulationPoisson::rhs(int equationI,
			       const GroupOfDof& god) const{
    
  // Init //
  fullMatrix<double>  jac(3, 3);        
  MElement& element = const_cast<MElement&>(god.getGeoElement());
  double integral   = 0;

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = element.getJacobian((*gC)(g, 0), 
				     (*gC)(g, 1), 
				     (*gC)(g, 2), 
				     jac);

    double phi = (*basis)[equationI].at((*gC)(g, 0), 
					(*gC)(g, 1), 
					(*gC)(g, 2));

    integral += 1 * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}
