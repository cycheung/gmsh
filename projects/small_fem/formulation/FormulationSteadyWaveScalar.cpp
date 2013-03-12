#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "GroupOfJacobian.h"

#include "Exception.h"
#include "FormulationSteadyWaveScalar.h"

using namespace std;

// Pi  = atan(1) * 4
// Mu  = 4 * Pi * 10^-7
// Eps = 8.85418781762 * 10^−12
//const double FormulationSteadyWaveScalar::mu  = 4 * atan(1) * 4 * 1E-7;
//const double FormulationSteadyWaveScalar::eps = 8.85418781762E-12;

const double FormulationSteadyWaveScalar::mu  = 1;
const double FormulationSteadyWaveScalar::eps = 1;

FormulationSteadyWaveScalar::FormulationSteadyWaveScalar(GroupOfElement& goe,
							 double k,
							 unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Scalar SteadyWave formulation of order 0");

  // Wave Number Squared //
  kSquare = k * k;

  // Function Space & Basis//
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature Data (Term One) //
  // NB: We need to integrad a grad * grad !
  //     and order(rot f) = order(f) - 1
  fullMatrix<double> gC1;
  fullVector<double> gW1;

  // Gaussian Quadrature Data (Term Two) //
  // NB: We need to integrad a f * f !
  fullMatrix<double> gC2;
  fullVector<double> gW2;

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), 2 * (order - 1), gC1, gW1);
  gaussIntegration::get(goe.get(0).getType(), 2 *  order     , gC2, gW2);

  // Local Terms //
  basis->preEvaluateDerivatives(gC1);
  basis->preEvaluateFunctions(gC2);

  GroupOfJacobian jac1(goe, gC1, "invert");
  GroupOfJacobian jac2(goe, gC2, "jacobian");

  localTerms1 = new TermGradGrad(jac1, *basis, gW1);
  localTerms2 = new TermFieldField(jac2, *basis, gW2);
}

FormulationSteadyWaveScalar::~FormulationSteadyWaveScalar(void){
  delete basis;
  delete fspace;

  delete localTerms1;
  delete localTerms2;
}

double FormulationSteadyWaveScalar::weak(unsigned int dofI, unsigned int dofJ,
                                         unsigned int elementId) const{
  return
    localTerms1->getTerm(dofI, dofJ, elementId) / mu -
    localTerms2->getTerm(dofI, dofJ, elementId) * eps * kSquare;
}

double FormulationSteadyWaveScalar::rhs(unsigned int equationI,
                                        unsigned int elementId) const{
  return 0;
}

bool FormulationSteadyWaveScalar::isGeneral(void) const{
  return false;
}

double FormulationSteadyWaveScalar::weakB(unsigned int dofI,
                                          unsigned int dofJ,
                                          unsigned int elementId) const{
  return 0;
}

const FunctionSpace& FormulationSteadyWaveScalar::fs(void) const{
  return *fspace;
}
