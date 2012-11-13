#include <cmath>

#include "GaussIntegration.h"
#include "Mapper.h"
#include "Polynomial.h"

#include "FormulationProjectionVector.h"

using namespace std;

FormulationProjectionVector::FormulationProjectionVector(const GroupOfElement& goe,
							 fullVector<double> (*f)(fullVector<double>& xyz),
							 unsigned int order){
  // Vector to Project //
  this->f = f;

  // Gaussian Quadrature Data  // 
  // NB: We need to integrad f_i * f_j or f_i * g
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)

  // if Order == 0 --> we want Nedelec Basis of ordre *almost* one //
  if(order != 0)
    gaussIntegration::get(goe.get(0).getType(), order + order, *gC, *gW);

  else
    gaussIntegration::get(goe.get(0).getType(), 1 + 1, *gC, *gW);
  
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
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double> invJac(3, 3);       

  double integral = 0;
  double det;

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
  fullVector<double> phi(3);
  double det;

  fullVector<double> xyz(3);
  SPoint3            pxyz;
  fullVector<double> fxyz;

  double integral = 0;
  fullMatrix<double> invJac(3, 3);       

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const vector<const vector<Polynomial>*> fun = 
    fspace->getLocalFunctions(element);  

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){  
    // Compute phi 
    det = celement.getJacobian((*gC)(g, 0), 
			       (*gC)(g, 1), 
			       (*gC)(g, 2), 
			       invJac);
    invJac.invertInPlace();
  

    phi = Mapper::grad(Polynomial::at(*fun[equationI],
				      (*gC)(g, 0), 
				      (*gC)(g, 1),
				      (*gC)(g, 2)),
		       invJac);
 

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
