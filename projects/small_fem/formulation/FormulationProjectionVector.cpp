#include <cmath>

#include "GaussIntegration.h"
#include "Mapper.h"
#include "Polynomial.h"

#include "FormulationProjectionVector.h"

using namespace std;

FormulationProjectionVector::
FormulationProjectionVector(fullVector<double> (*f)(fullVector<double>& xyz),
			    FunctionSpaceEdge& fs){
  // Vector to Project //
  this->f = f;

  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(0);

  // Gaussian Quadrature Data  //
  // NB: We need to integrad f_i * f_j or f_i * g
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(fs.getSupport().get(0).getType(), 2 * basis->getOrder(), *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // PreEvaluate
  basis->preEvaluateFunctions(*gC);
}

FormulationProjectionVector::~FormulationProjectionVector(void){
  delete gC;
  delete gW;
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

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    det = celement.getJacobian((*gC)(g, 0),
			       (*gC)(g, 1),
			       (*gC)(g, 2),
			       invJac);
    invJac.invertInPlace();

    phiI = Mapper::grad(eFun(dofI, g * 3),
			eFun(dofI, g * 3 + 1),
			eFun(dofI, g * 3 + 2),
			invJac);

    phiJ = Mapper::grad(eFun(dofJ, g * 3),
			eFun(dofJ, g * 3 + 1),
			eFun(dofJ, g * 3 + 2),
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

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    // Compute phi
    det = celement.getJacobian((*gC)(g, 0),
			       (*gC)(g, 1),
			       (*gC)(g, 2),
			       invJac);
    invJac.invertInPlace();

    phi = Mapper::grad(eFun(equationI, g * 3),
		       eFun(equationI, g * 3 + 1),
		       eFun(equationI, g * 3 + 2),
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
