#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "GroupOfElement.h"
#include "Quadrature.h"

#include "FormulationProjectionVector.h"

using namespace std;

FormulationProjectionVector::
FormulationProjectionVector(fullVector<double> (*f)(fullVector<double>& xyz),
			    FunctionSpaceVector& fs){
  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(0);

  // Domain //
  GroupOfElement& goe = fs.getSupport();

  // Gaussian Quadrature //
  Quadrature gauss(goe.get(0).getType(), basis->getOrder(), 2);

  const fullMatrix<double>& gC = gauss.getPoints();
  const fullVector<double>& gW = gauss.getWeights();

  // Local Terms //
  basis->preEvaluateFunctions(gC);

  GroupOfJacobian jac(goe, gC, "invert");

  localTerms1 = new TermGradGrad(jac, *basis, gW);
  localTerms2 = new TermProjectionGrad(jac, *basis, gW, gC, f);
}

FormulationProjectionVector::~FormulationProjectionVector(void){
  delete localTerms1;
  delete localTerms2;
}

double FormulationProjectionVector::weak(unsigned int dofI, unsigned int dofJ,
                                         unsigned int elementId) const{

  return localTerms1->getTerm(dofI, dofJ, elementId);
}

double FormulationProjectionVector::rhs(unsigned int equationI,
					unsigned int elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}

bool FormulationProjectionVector::isGeneral(void) const{
  return false;
}

double FormulationProjectionVector::weakB(unsigned int dofI,
                                          unsigned int dofJ,
                                          unsigned int elementId) const{
  return 0;
}

const FunctionSpace& FormulationProjectionVector::fs(void) const{
  return *fspace;
}
