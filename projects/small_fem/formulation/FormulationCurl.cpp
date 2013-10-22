#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationCurl.h"

using namespace std;

FormulationCurl::FormulationCurl(GroupOfElement& goe,
                                 size_t order){
  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    1, order, "hierarchical");

  fspace = new FunctionSpaceVector(goe, *basis);

  // Gaussian Quadrature //
  Quadrature gauss(goe.get(0).getType(), order - 1, 2);

  const fullMatrix<double>& gC = gauss.getPoints();
  const fullVector<double>& gW = gauss.getWeights();

  // Local Terms //
  basis->preEvaluateDerivatives(gC);

  GroupOfJacobian jac(goe, *basis, gC, "jacobian");

  localTerms = new TermCurlCurl(jac, *basis, gW);
}

FormulationCurl::~FormulationCurl(void){
  delete basis;
  delete fspace;
  delete localTerms;
}

double FormulationCurl::weak(size_t dofI, size_t dofJ,
                                size_t elementId) const{

  return localTerms->getTerm(dofI, dofJ, elementId);
}


double FormulationCurl::rhs(size_t equationI,
                               size_t elementId) const{
  return 0;
}

bool FormulationCurl::isGeneral(void) const{
  return false;
}

double FormulationCurl::weakB(size_t dofI,
                                 size_t dofJ,
                                 size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationCurl::fs(void) const{
  return *fspace;
}
