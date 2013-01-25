#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Jacobian.h"
#include "GroupOfElement.h"

#include "FormulationProjectionVector.h"

using namespace std;

FormulationProjectionVector::
FormulationProjectionVector(fullVector<double> (*f)(fullVector<double>& xyz),
			    FunctionSpaceVector& fs){
  // Save f //
  this->f = f;

  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(0);

  // Domain //
  GroupOfElement& goe = fs.getSupport();

  // Gaussian Quadrature Data  //
  // NB: We need to integrad f_i * f_j or f_i * g
  fullMatrix<double> gC;
  fullVector<double> gW;

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), 2 * basis->getOrder(), gC, gW);

  // Local Terms //
  basis->preEvaluateFunctions(gC);
  goe.orientAllElements(*basis);

  Jacobian jac(goe, gC);
  jac.computeInvertJacobians();

  localTerms1 = new TermHCurl(jac, *basis, gW);
  localTerms2 = new TermProjectionHCurl(jac, *basis, gW, gC, f);
}

FormulationProjectionVector::~FormulationProjectionVector(void){
  delete localTerms1;
  delete localTerms2;
}

double FormulationProjectionVector::weak(unsigned int dofI, unsigned int dofJ,
                                         const GroupOfDof& god) const{

  return localTerms1->getTerm(dofI, dofJ, god);
}

double FormulationProjectionVector::rhs(unsigned int equationI,
					const GroupOfDof& god) const{

  return localTerms2->getTerm(0, equationI, god);
}

const FunctionSpace& FormulationProjectionVector::fs(void) const{
  return *fspace;
}
