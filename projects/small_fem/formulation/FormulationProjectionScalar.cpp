#include <cmath>

#include "GaussIntegration.h"
#include "Mapper.h"
#include "Polynomial.h"

#include "Exception.h"

#include "FormulationProjectionScalar.h"

using namespace std;

FormulationProjectionScalar::
FormulationProjectionScalar(double (*f)(fullVector<double>& xyz),
			    const FunctionSpaceNode& fs){
  // Save f //
  this->f = f;

  // Save fspace //
  fspace = &fs;

  // Gaussian Quadrature Data  // 
  // NB: We need to integrad f_i * f_j or f_i * g
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(fs.getSupport().get(0).getType(), 2 * fs.getOrder(), *gC, *gW);

  G = gW->size(); // Nbr of Gauss points
}

FormulationProjectionScalar::~FormulationProjectionScalar(void){
  delete gC;
  delete gW;
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
  double phi;
  double det;

  fullVector<double> xyz(3);
  SPoint3            pxyz;
  double             fxyz;

  double integral = 0;
  fullMatrix<double> jac(3, 3);

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const Polynomial*> fun = 
    fspace->getLocalFunctions(element);  

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    // Compute phi 
    det = celement.getJacobian((*gC)(g, 0), 
			       (*gC)(g, 1), 
			       (*gC)(g, 2), 
			       jac);

    phi = fun[equationI]->at((*gC)(g, 0), 
			     (*gC)(g, 1),
			     (*gC)(g, 2));
    
    // Compute f in the *physical* coordinate
    celement.pnt((*gC)(g, 0), 
		 (*gC)(g, 1), 
		 (*gC)(g, 2), 
		 pxyz);
    
    xyz(0) = pxyz.x();
    xyz(1) = pxyz.y();
    xyz(2) = pxyz.z();
    
    fxyz = f(xyz);

    // Integrate
    integral += fxyz * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}
