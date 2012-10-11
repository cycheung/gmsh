#include <cmath>

#include "GaussIntegration.h"
#include "Mapper.h"
#include "Polynomial.h"

#include "MVertex.h"
#include "Exception.h"

#include "FormulationProjectionVector.h"

using namespace std;

FormulationProjectionVector::FormulationProjectionVector(const GroupOfElement& goe,
							 fullVector<double> (*f)(fullVector<double>& xyz),
							 unsigned int order){
  // Vector to Project //
  this->f = f;

  // Can't have 0th order //
  if(order == 0)
    throw 
      Exception("Can't have a Projection of order 0");

  // Gaussian Quadrature Data  // 
  // NB: We need to integrad f_i * f_j or f_i * g
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), order + order, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Function Space //
  fspace = new FunctionSpaceEdge(goe, order);
}

FormulationProjectionVector::~FormulationProjectionVector(void){
  delete gC;
  delete gW;
  delete fspace;
}

double FormulationProjectionVector::weak(int dofI, int dofJ, 
					 const GroupOfDof& god) const{
  // Init //
  double det;
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
    det = celement.getJacobian((*gC)(g, 0), 
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

double FormulationProjectionVector::rhs(int equationI,
					const GroupOfDof& god) const{
  // Init //
  fullMatrix<double>  invJac(3, 3);        
  fullMatrix<double>     jac(3, 3);        
  fullVector<double>  phi(3);
  double integral = 0;
  double det;

  MVertex* vertex;
  fullVector<double>  uvw(3);
  fullVector<double> oxyz(3);
  fullVector<double>  xyz(3);
  fullVector<double> fxyz(3);


  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const vector<Polynomial>*> fun = 
    fspace->getLocalFunctions(element);  

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){  
    det = celement.getJacobian((*gC)(g, 0), 
			       (*gC)(g, 1), 
			       (*gC)(g, 2), 
			       jac);
    jac.invert(invJac);

    // Parametric coordinate 
    uvw(0) = (*gC)(g, 0);
    uvw(1) = (*gC)(g, 1);
    uvw(2) = (*gC)(g, 2);
  
    // Compute phi 
    phi = Mapper::grad(Polynomial::at(*fun[equationI],
				      uvw(0), 
				      uvw(1),
				      uvw(2)),
		       invJac);
 
    // Compute f in the *physical* coordinate
    //  --> Get *physical* coordinate
    //       --> Get Origin Point of Element
    vertex = celement.getVertex(0);
    oxyz(0) = vertex->x();
    oxyz(1) = vertex->y();
    oxyz(2) = vertex->z();
    
    //       --> Map
    xyz = Mapper::map(uvw, oxyz, jac);

    // --> Evaluate f
    fxyz = f(xyz);

    // Interate
    integral += fxyz * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}
