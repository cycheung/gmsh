#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "GroupOfElement.h"
#include "Quadrature.h"

#include "FormulationProjectionScalar.h"

using namespace std;

FormulationProjectionScalar::
FormulationProjectionScalar(double (*f)(fullVector<double>& xyz),
                            FunctionSpaceScalar& fs){
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
  GroupOfJacobian jac(goe, *basis, gC, "jacobian");

  localTerms1 = new TermFieldField(jac, *basis, gW);
  localTerms2 = new TermProjectionField(jac, *basis, gW, gC, f);
}

FormulationProjectionScalar::~FormulationProjectionScalar(void){
  delete localTerms2;
  delete localTerms1;
}

double FormulationProjectionScalar::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{

  return localTerms1->getTerm(dofI, dofJ, elementId);
}

double FormulationProjectionScalar::rhs(size_t equationI,
                                        size_t elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}

bool FormulationProjectionScalar::isGeneral(void) const{
  return false;
}

double FormulationProjectionScalar::weakB(size_t dofI,
                                          size_t dofJ,
                                          size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationProjectionScalar::fs(void) const{
  return *fspace;
}
