#include <cmath>

#include "GaussIntegration.h"
#include "Mapper.h"
#include "Polynomial.h"

#include "MVertex.h"
#include "Exception.h"

#include "FormulationProjectionScalar.h"

using namespace std;

FormulationProjectionScalar::FormulationProjectionScalar(const GroupOfElement& goe,
							 double (*f)(fullVector<double>& xyz),
							 unsigned int ordre){
  // Save f //
  this->f = f;

  // Can't have 0th order //
  if(ordre == 0)
    throw 
      Exception("Can't have a Projection of order 0");

  // Gaussian Quadrature Data  // 
  // NB: We need to integrad f_i * f_j or f_i * g
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), ordre + ordre, *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // Function Space //
  fspace = new FunctionSpaceNode(goe, ordre);
}

FormulationProjectionScalar::~FormulationProjectionScalar(void){
  delete gC;
  delete gW;
  delete fspace;
}

double FormulationProjectionScalar::weak(int dofI, int dofJ, 
				   const GroupOfDof& god) const{
  // Init //
  double det;
  double phiI;
  double phiJ;
  fullMatrix<double> jac(3, 3);   
  double integral = 0;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const Polynomial*> fun = 
    fspace->getLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    det = celement.getJacobian((*gC)(g, 0), 
			       (*gC)(g, 1), 
			       (*gC)(g, 2), 
			       jac);

    phiI = fun[dofI]->at((*gC)(g, 0),
			 (*gC)(g, 1),
			 (*gC)(g, 2));
  
    phiJ = fun[dofJ]->at((*gC)(g, 0), 
			 (*gC)(g, 1),
			 (*gC)(g, 2));

    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return integral;
}

double FormulationProjectionScalar::rhs(int equationI,
					const GroupOfDof& god) const{
  // Init //
  double det;
  double phi;

  fullVector<double>  uvw(3);
  fullVector<double> oxyz(3);
  fullVector<double>  xyz(3);
  double fxyz;

  double integral = 0;
  fullMatrix<double> jac(3, 3);
  MVertex* vertex;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const Polynomial*> fun = 
    fspace->getLocalFunctions(element);  

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    // Parametric coordinate 
    uvw(0) = (*gC)(g, 0);
    uvw(1) = (*gC)(g, 1);
    uvw(2) = (*gC)(g, 2);
  
    // Compute phi 
    det = celement.getJacobian(uvw(0), 
			       uvw(1), 
			       uvw(2), 
			       jac);

    phi = fun[equationI]->at(uvw(0), 
			     uvw(1),
			     uvw(2));
    
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
