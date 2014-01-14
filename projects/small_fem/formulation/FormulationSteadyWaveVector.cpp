#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationSteadyWaveVector.h"

using namespace std;

FormulationSteadyWaveVector::FormulationSteadyWaveVector(GroupOfElement& goe,
                                                         double k,
                                                         size_t order){
  // Wavenumber Squared //
  kSquare = k * k;

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    1, order, "hierarchical");

  fspace = new FunctionSpaceVector(goe, *basis);

  // Gaussian Quadrature //
  Quadrature gaussCurlCurl(goe.get(0).getType(), order - 1, 2);
  Quadrature gaussFF(goe.get(0).getType(), order, 2);

  const fullMatrix<double>& gC1 = gaussCurlCurl.getPoints();
  const fullVector<double>& gW1 = gaussCurlCurl.getWeights();

  const fullMatrix<double>& gC2 = gaussFF.getPoints();
  const fullVector<double>& gW2 = gaussFF.getWeights();

  // Local Terms //
  basis->preEvaluateDerivatives(gC1);
  basis->preEvaluateFunctions(gC2);

  GroupOfJacobian jac1(goe, *basis, gC1, "jacobian");
  GroupOfJacobian jac2(goe, *basis, gC2, "invert");

  localTerms1 = new TermCurlCurl(jac1, *basis, gW1);
  localTerms2 = new TermGradGrad(jac2, *basis, gW2);
}

FormulationSteadyWaveVector::~FormulationSteadyWaveVector(void){
  delete basis;
  delete fspace;

  delete localTerms1;
  delete localTerms2;
}

double FormulationSteadyWaveVector::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{
  return
    localTerms1->getTerm(dofI, dofJ, elementId) -
    localTerms2->getTerm(dofI, dofJ, elementId) * kSquare;
}

double FormulationSteadyWaveVector::rhs(size_t equationI,
                                        size_t elementId) const{
  return 0;
}

bool FormulationSteadyWaveVector::isGeneral(void) const{
  return false;
}

double FormulationSteadyWaveVector::weakB(size_t dofI,
                                          size_t dofJ,
                                          size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationSteadyWaveVector::fs(void) const{
  return *fspace;
}
