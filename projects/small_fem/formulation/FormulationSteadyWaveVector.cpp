#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Timer.h"
#include "FormulationSteadyWaveVector.h"

using namespace std;

// Pi  = atan(1) * 4
// Mu  = 4 * Pi * 10^-7
// Eps = 8.85418781762 * 10^−12
//const double FormulationSteadyWaveVector::mu  = 4 * atan(1) * 4 * 1E-7;
//const double FormulationSteadyWaveVector::eps = 8.85418781762E-12;

const double FormulationSteadyWaveVector::mu  = 1;
const double FormulationSteadyWaveVector::eps = 1;

FormulationSteadyWaveVector::FormulationSteadyWaveVector(GroupOfElement& goe,
                                                         double k,
                                                         unsigned int order){
  Timer timer, timerGoj;
  timer.start();

  // Wave Number Squared //
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

  timerGoj.start();
  GroupOfJacobian jac1(goe, *basis, gC1, "jacobian");
  GroupOfJacobian jac2(goe, *basis, gC2, "invert");
  timerGoj.stop();


  localTerms1 = new TermCurlCurl(jac1, *basis, gW1);
  localTerms2 = new TermGradGrad(jac2, *basis, gW2);

  timer.stop();

  cout << "Jacs: " << timerGoj.time() << " " << timerGoj.unit() << endl;
  cout << "Full: " << timer.time() << " " << timer.unit() << endl;
}

FormulationSteadyWaveVector::~FormulationSteadyWaveVector(void){
  delete basis;
  delete fspace;

  delete localTerms1;
  delete localTerms2;
}

double FormulationSteadyWaveVector::weak(unsigned int dofI, unsigned int dofJ,
                                         unsigned int elementId) const{
  return
    localTerms1->getTerm(dofI, dofJ, elementId) / mu -
    localTerms2->getTerm(dofI, dofJ, elementId) * eps * kSquare;
}

double FormulationSteadyWaveVector::rhs(unsigned int equationI,
                                        unsigned int elementId) const{
  return 0;
}

bool FormulationSteadyWaveVector::isGeneral(void) const{
  return false;
}

double FormulationSteadyWaveVector::weakB(unsigned int dofI,
                                          unsigned int dofJ,
                                          unsigned int elementId) const{
  return 0;
}

const FunctionSpace& FormulationSteadyWaveVector::fs(void) const{
  return *fspace;
}
