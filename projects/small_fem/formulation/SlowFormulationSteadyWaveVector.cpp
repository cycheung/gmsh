#include <cmath>

#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Mapper.h"

#include "Timer.h"
#include "SlowFormulationSteadyWaveVector.h"

using namespace std;

// Pi  = atan(1) * 4
// Mu  = 4 * Pi * 10^-7
// Eps = 8.85418781762 * 10^−12
//const double SlowFormulationSteadyWaveVector::mu  = 4 * atan(1) * 4 * 1E-7;
//const double SlowFormulationSteadyWaveVector::eps = 8.85418781762E-12;

const double SlowFormulationSteadyWaveVector::mu  = 1;
const double SlowFormulationSteadyWaveVector::eps = 1;

SlowFormulationSteadyWaveVector::SlowFormulationSteadyWaveVector(GroupOfElement& goe,
                                                                 double k,
                                                                 unsigned int order){
  // Wave Number Squared //
  kSquare = k * k;

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    1, order, "hierarchical");

  fspace = new FunctionSpaceVector(goe, *basis);

  // Gaussian Quadrature Data (Term One) //
  // NB: We need to integrad a rot * rot !
  //     and order(rot f) = order(f) - 1
  gC1 = new fullMatrix<double>();
  gW1 = new fullVector<double>();

  // Gaussian Quadrature Data (Term Two) //
  // NB: We need to integrad a f * f !
  gC2 = new fullMatrix<double>();
  gW2 = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)

  // if Order == 0 --> we want Nedelec Basis of ordre *almost* one //
  if(order == 0){
    gaussIntegration::get(goe.get(0).getType(), 0, *gC1, *gW1);
    gaussIntegration::get(goe.get(0).getType(), 2, *gC2, *gW2);
  }

  else{
    gaussIntegration::get(goe.get(0).getType(), (order - 1) + (order - 1), *gC1, *gW1);
    gaussIntegration::get(goe.get(0).getType(), order + order, *gC2, *gW2);
  }

  // Nbr of Gauss points
  G1 = gW1->size();
  G2 = gW2->size();

  // PreEvaluate Functions //
  Timer timer;
  timer.start();

  basis->preEvaluateDerivatives(*gC1);
  basis->preEvaluateFunctions(*gC2);

  // PreEvaluate Jacobians //
  goe.orientAllElements(*basis);

  jac1 = new Jacobian(goe, *gC1);
  jac2 = new Jacobian(goe, *gC2);

  jac1->computeJacobians();
  jac2->computeInvertJacobians();

  timer.stop();
  cout << "Precomputing Time: " << timer.time() << " " << timer.unit() << endl;
}

SlowFormulationSteadyWaveVector::~SlowFormulationSteadyWaveVector(void){
  delete gC1;
  delete gW1;
  delete gC2;
  delete gW2;
  delete basis;
  delete fspace;
  delete jac1;
  delete jac2;
}

double SlowFormulationSteadyWaveVector::weak(unsigned int dofI,unsigned int dofJ,
                                             const GroupOfDof& god) const{
  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);

  const fullMatrix<double>* jac;

  double integral1 = 0;
  double integral2 = 0;
  double det;

  // Get Element //
  const MElement& element = god.getGeoElement();

  // Get Basis Functions //
  const fullMatrix<double>& eCurlFun =
    basis->getPreEvaluatedDerivatives(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& MJac =
    jac1->getJacobian(element);

  const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
    jac2->getInvertJacobian(element);

  // Loop over Integration Point (Term 1) //
  for(int g = 0; g < G1; g++){
    det = MJac[g]->second;
    jac = MJac[g]->first;

    Mapper::hDiv(eCurlFun, dofI, g, *jac, det, phiI);
    Mapper::hDiv(eCurlFun, dofJ, g, *jac, det, phiJ);

    integral1 +=
      ((phiI * phiJ) / mu) * fabs(det) * (*gW1)(g);
  }


  // Loop over Integration Point (Term 2) //
  for(int g = 0; g < G2; g++){
    det = invJac[g]->second;
    jac = invJac[g]->first;

    Mapper::hCurl(eFun, dofI, g, *jac, phiI);
    Mapper::hCurl(eFun, dofJ, g, *jac, phiJ);

    integral2 +=
      ((phiI * phiJ) * eps * kSquare) * fabs(det) * (*gW2)(g);
  }

  return integral1 - integral2;
}

double SlowFormulationSteadyWaveVector::rhs(unsigned int equationI,
                                            const GroupOfDof& god) const{
  return 0;
}

const FunctionSpace& SlowFormulationSteadyWaveVector::fs(void) const{
  return *fspace;
}

