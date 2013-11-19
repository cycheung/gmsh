#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationSteadyWaveScalar.h"

using namespace std;

const double FormulationSteadyWaveScalar::cSquare = 1;

FormulationSteadyWaveScalar::FormulationSteadyWaveScalar(GroupOfElement& goe,
                                                         double omega,
                                                         size_t order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Scalar SteadyWave formulation of order 0");

  // Pulsation Squared //
  omegaSquare = omega * omega;

  // Function Space & Basis//
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature //
  Quadrature gaussGradGrad(goe.get(0).getType(), order - 1, 2);
  Quadrature gaussFF(goe.get(0).getType(), order, 2);

  const fullMatrix<double>& gC1 = gaussGradGrad.getPoints();
  const fullVector<double>& gW1 = gaussGradGrad.getWeights();

  const fullMatrix<double>& gC2 = gaussFF.getPoints();
  const fullVector<double>& gW2 = gaussFF.getWeights();

  // Local Terms //
  basis->preEvaluateDerivatives(gC1);
  basis->preEvaluateFunctions(gC2);

  GroupOfJacobian jac1(goe, *basis, gC1, "invert");
  GroupOfJacobian jac2(goe, *basis, gC2, "jacobian");

  localTerms1 = new TermGradGrad(jac1, *basis, gW1);
  localTerms2 = new TermFieldField(jac2, *basis, gW2);
}

FormulationSteadyWaveScalar::~FormulationSteadyWaveScalar(void){
  delete basis;
  delete fspace;

  delete localTerms1;
  delete localTerms2;
}

double FormulationSteadyWaveScalar::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{
  return
    localTerms1->getTerm(dofI, dofJ, elementId) -
    localTerms2->getTerm(dofI, dofJ, elementId) * omegaSquare / cSquare;
}

double FormulationSteadyWaveScalar::rhs(size_t equationI,
                                        size_t elementId) const{
  return 0;
}

bool FormulationSteadyWaveScalar::isGeneral(void) const{
  return false;
}

double FormulationSteadyWaveScalar::weakB(size_t dofI,
                                          size_t dofJ,
                                          size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationSteadyWaveScalar::fs(void) const{
  return *fspace;
}
