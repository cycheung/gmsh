#include <cmath>

#include "Exception.h"
#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "BasisScalar.h"
#include "Mapper.h"

#include "FunctionSpaceNode.h"
#include "FormulationPoisson.h"

//#include <iostream>

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
     basis = &(base.getFunctions(0));
  revBasis = &(base.getFunctions(1));

  // Take gradient
  unsigned int basisSize = basis->size();

     gradBasis = new vector<Polynomial>[basisSize];
  revGradBasis = new vector<Polynomial>[basisSize];

  for(unsigned int i = 0; i < basisSize; i++){
       gradBasis[i] = (*basis)[i]->gradient();
       revGradBasis[i] = (*revBasis)[i]->gradient();
  }
}

FormulationPoisson::~FormulationPoisson(void){
  delete   gC;
  delete   gW;
  delete   fspace;
  delete[] gradBasis;
  delete[] revGradBasis;
}

double FormulationPoisson::weak(int nodeI, int nodeJ, 
				const GroupOfDof& god) const{

  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double>  invJac(3, 3);        
  MElement& element = const_cast<MElement&>(god.getGeoElement());
  double integral   = 0;

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = element.getJacobian((*gC)(g, 0), 
				     (*gC)(g, 1), 
				     (*gC)(g, 2), 
				     invJac);
    /*    
    cout << "Element: " << element.getNum() << endl;
    cout << "Origin : [" 
	 << element.getVertex(0)->x() << ", "
	 << element.getVertex(0)->y() << ", "
	 << element.getVertex(0)->z() << "]" << endl;

    invJac.print();
    */
    invJac.invertInPlace();

    if(god.getOrientation(nodeI) == - 1)
      phiI = Mapper::grad(Polynomial::at(revGradBasis[nodeI], 
					 (*gC)(g, 0), 
					 (*gC)(g, 1),
					 (*gC)(g, 2)),
			  invJac);
    
    else
      phiI = Mapper::grad(Polynomial::at(gradBasis[nodeI], 
					 (*gC)(g, 0), 
					 (*gC)(g, 1),
					 (*gC)(g, 2)),
			  invJac);

    
    if(god.getOrientation(nodeJ) == - 1)
      phiJ = Mapper::grad(Polynomial::at(revGradBasis[nodeJ], 
					 (*gC)(g, 0), 
					 (*gC)(g, 1), 
					 (*gC)(g, 2)),
			  invJac);
	
    else
      phiJ = Mapper::grad(Polynomial::at(gradBasis[nodeJ], 
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
  double phi;

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = element.getJacobian((*gC)(g, 0), 
				     (*gC)(g, 1), 
				     (*gC)(g, 2), 
				     jac);

    if(god.getOrientation(equationI) == -1)
      phi = (*revBasis)[equationI]->at((*gC)(g, 0), 
				       (*gC)(g, 1), 
				       (*gC)(g, 2));

    else
      phi = (*basis)[equationI]->at((*gC)(g, 0), 
				    (*gC)(g, 1), 
				    (*gC)(g, 2));    
    

    integral += 1 * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}
