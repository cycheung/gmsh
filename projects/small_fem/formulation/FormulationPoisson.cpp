#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationPoisson.h"

using namespace std;

// Poisson //
FormulationPoisson::FormulationPoisson(GroupOfElement& goe,
                                       double (*f)(fullVector<double>& xyz),
                                       size_t order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Poisson formulation of order 0");

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Source Term //
  fSource = f;

  // Gaussian Quadrature //
  Quadrature gaussGradGrad(goe.get(0).getType(), order - 1, 2);
  Quadrature gaussFF(goe.get(0).getType(), order, 2);

  const fullMatrix<double>& gCL = gaussGradGrad.getPoints();
  const fullVector<double>& gWL = gaussGradGrad.getWeights();

  const fullMatrix<double>& gCR = gaussFF.getPoints();
  const fullVector<double>& gWR = gaussFF.getWeights();

  // Local Terms //
  basis->preEvaluateDerivatives(gCL);
  basis->preEvaluateFunctions(gCR);

  GroupOfJacobian jacL(goe, *basis, gCL, "invert");
  GroupOfJacobian jacR(goe, *basis, gCR, "jacobian");

  localTermsL = new TermGradGrad(jacL, *basis, gWL);
  localTermsR = new TermProjectionField(jacR, *basis, gWR, gCR, fSource);
}

FormulationPoisson::~FormulationPoisson(void){
  delete basis;
  delete fspace;

  delete localTermsL;
  delete localTermsR;
}

double FormulationPoisson::weak(size_t dofI, size_t dofJ,
                                size_t elementId) const{

  return localTermsL->getTerm(dofI, dofJ, elementId);
}

double FormulationPoisson::rhs(size_t equationI,
                               size_t elementId) const{

  return localTermsR->getTerm(0, equationI, elementId);
}

bool FormulationPoisson::isGeneral(void) const{
  return false;
}

double FormulationPoisson::weakB(size_t dofI,
                                 size_t dofJ,
                                 size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationPoisson::fs(void) const{
  return *fspace;
}
