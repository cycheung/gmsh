#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationLaplace.h"

using namespace std;

FormulationLaplace::FormulationLaplace(GroupOfElement& goe,
                                       size_t order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Laplace formulation of order 0");

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature //
  Quadrature gauss(goe.get(0).getType(), order - 1, 2);

  const fullMatrix<double>& gC = gauss.getPoints();
  const fullVector<double>& gW = gauss.getWeights();

  // Local Terms //
  basis->preEvaluateDerivatives(gC);

  GroupOfJacobian jac(goe, *basis, gC, "invert");

  localTerms = new TermGradGrad(jac, *basis, gW);
}

FormulationLaplace::~FormulationLaplace(void){
  delete basis;
  delete fspace;
  delete localTerms;
}

double FormulationLaplace::weak(size_t dofI, size_t dofJ,
                                size_t elementId) const{

  return localTerms->getTerm(dofI, dofJ, elementId);
}


double FormulationLaplace::rhs(size_t equationI,
                               size_t elementId) const{
  return 0;
}

bool FormulationLaplace::isGeneral(void) const{
  return false;
}

double FormulationLaplace::weakB(size_t dofI,
                                 size_t dofJ,
                                 size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationLaplace::fs(void) const{
  return *fspace;
}
