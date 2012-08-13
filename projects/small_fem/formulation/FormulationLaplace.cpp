#include <cmath>

#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "BasisScalar.h"
#include "Mapper.h"

#include "FunctionSpaceNode.h"
#include "FormulationLaplace.h"

using namespace std;

FormulationLaplace::FormulationLaplace(const GroupOfElement& goe){
  // Gaussian Quadrature Data //
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), 2, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Function Space //
  FunctionSpaceNode* fspace = new FunctionSpaceNode(goe, 1);
  this->fspace              = fspace;

  // Basis //
  // Get Basis
  const BasisScalar& base = fspace->getBasis(goe.get(0));
  const vector<Polynomial>& basis = base.getFunctions();

  // Take gradient
  unsigned int basisSize = basis.size();

  gradBasis = new vector<Polynomial>[basisSize];

  for(unsigned int i = 0; i < basisSize; i++)
    gradBasis[i] = basis[i].gradient();  
}

FormulationLaplace::~FormulationLaplace(void){
  delete   gC;
  delete   gW;
  delete   fspace;
  delete[] gradBasis;
}
//#include <cstdio>
double FormulationLaplace::weak(int nodeI, int nodeJ, 
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
    /*
    printf("%d\n", element.getNum());
    printf("zero : (%lf\t%lf\t%lf)\n", element.getVertex(0)->x(),
	   element.getVertex(0)->y(),
	   element.getVertex(0)->z());
    printf("one  : (%lf\t%lf\t%lf)\n", element.getVertex(1)->x(),
	   element.getVertex(1)->y(),
	   element.getVertex(1)->z());
    printf("two  : (%lf\t%lf\t%lf)\n", element.getVertex(2)->x(),
	   element.getVertex(2)->y(),
	   element.getVertex(2)->z());
    printf("three: (%lf\t%lf\t%lf)\n", element.getVertex(3)->x(),
	   element.getVertex(3)->y(),
	   element.getVertex(3)->z());
 
    invJac.print();
    
    printf("(%lf\t%lf\t%lf)\n", (*gC)(g, 0), (*gC)(g, 1), (*gC)(g, 2));
    */
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
