#include <cmath>

#include "BasisGenerator.h"
#include "Quadrature.h"
#include "Mapper.h"

#include "FormulationSteadyWaveVectorSlow.h"

using namespace std;

const double FormulationSteadyWaveVectorSlow::cSquare = 1;

FormulationSteadyWaveVectorSlow::
FormulationSteadyWaveVectorSlow(GroupOfElement& goe,
                                double omega,
                                size_t order){
  // Pulsation Squared //
  omegaSquare = omega * omega;

  // Domain //
  this->goe = &goe;

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    1, order, "hierarchical");

  fspace = new FunctionSpaceVector(goe, *basis);

  // Gaussian Quadrature //
  Quadrature gaussCurlCurl(goe.get(0).getType(), order - 1, 2);
  Quadrature gaussFF(goe.get(0).getType(), order, 2);

  gC1 = new fullMatrix<double>(gaussCurlCurl.getPoints());
  gW1 = new fullVector<double>(gaussCurlCurl.getWeights());

  gC2 = new fullMatrix<double>(gaussFF.getPoints());
  gW2 = new fullVector<double>(gaussFF.getWeights());

  G1 = gW1->size();
  G2 = gW2->size();

  // PreEvaluate
  basis->preEvaluateDerivatives(*gC1);
  basis->preEvaluateFunctions(*gC2);

  jac1 = new GroupOfJacobian(goe, *basis, *gC1, "jacobian");
  jac2 = new GroupOfJacobian(goe, *basis, *gC2, "invert");
}

FormulationSteadyWaveVectorSlow::~FormulationSteadyWaveVectorSlow(void){
  delete gC1;
  delete gW1;
  delete gC2;
  delete gW2;
  delete jac1;
  delete jac2;
  delete basis;
  delete fspace;
}

double FormulationSteadyWaveVectorSlow::weak(size_t dofI,
                                             size_t dofJ,
                                             size_t elementId) const{
  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);

  const fullMatrix<double>* jac;

  double integral1 = 0;
  double integral2 = 0;
  double det;

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eCurlFun =
    basis->getPreEvaluatedDerivatives(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& MJac =
    jac1->getJacobian(elementId).getJacobianMatrix();

  const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
    jac2->getJacobian(elementId).getInvertJacobianMatrix();

  // Loop over Integration Point (Term 1) //
  for(int g = 0; g < G1; g++){
    det = MJac[g]->second;
    jac = MJac[g]->first;

    Mapper::hDiv(eCurlFun, dofI, g, *jac, det, phiI);
    Mapper::hDiv(eCurlFun, dofJ, g, *jac, det, phiJ);

    integral1 +=
      ((phiI * phiJ)) * fabs(det) * (*gW1)(g);
  }


  // Loop over Integration Point (Term 2) //
  for(int g = 0; g < G2; g++){
    det = invJac[g]->second;
    jac = invJac[g]->first;

    Mapper::hCurl(eFun, dofI, g, *jac, phiI);
    Mapper::hCurl(eFun, dofJ, g, *jac, phiJ);

    integral2 +=
      ((phiI * phiJ) * omegaSquare / cSquare) * fabs(det) * (*gW2)(g);
  }

  return integral1 - integral2;
}

double FormulationSteadyWaveVectorSlow::rhs(size_t equationI,
                                            size_t elementId) const{
  return 0;
}

bool FormulationSteadyWaveVectorSlow::isGeneral(void) const{
  return false;
}

double FormulationSteadyWaveVectorSlow::weakB(size_t dofI,
                                              size_t dofJ,
                                              size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationSteadyWaveVectorSlow::fs(void) const{
  return *fspace;
}
