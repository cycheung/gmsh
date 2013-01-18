#include <cmath>

#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Mapper.h"

#include "FormulationSteadyWaveVector.h"

using namespace std;

// Pi  = atan(1) * 4
// Mu  = 4 * Pi * 10^-7
// Eps = 8.85418781762 * 10^−12
//const double FormulationSteadyWaveVector::mu  = 4 * atan(1) * 4 * 1E-7;
//const double FormulationSteadyWaveVector::eps = 8.85418781762E-12;

const double FormulationSteadyWaveVector::mu  = 1;
const double FormulationSteadyWaveVector::eps = 1;

FormulationSteadyWaveVector::FormulationSteadyWaveVector(const GroupOfElement& goe,
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

  // PreEvaluate
  basis->preEvaluateDerivatives(*gC1);
  basis->preEvaluateFunctions(*gC2);
}

FormulationSteadyWaveVector::~FormulationSteadyWaveVector(void){
  delete gC1;
  delete gW1;
  delete gC2;
  delete gW2;
  delete basis;
  delete fspace;
}

double FormulationSteadyWaveVector::weak(int dofI, int dofJ,
					 const GroupOfDof& god) const{
  // Init Some Stuff //
  fullVector<double> curlPhiI(3);
  fullVector<double> curlPhiJ(3);
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);

  fullMatrix<double> jac(3, 3);
  fullMatrix<double> invJac(3, 3);

  double integral1 = 0;
  double integral2 = 0;
  double det;

  // Get Element and Basis Functions (+ Curl) //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);

  const fullMatrix<double>& eCurlFun =
    basis->getPreEvaluatedDerivatives(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Loop over Integration Point (Term 1) //
  for(int g = 0; g < G1; g++){
    det = celement.getJacobian((*gC1)(g, 0),
			       (*gC1)(g, 1),
			       (*gC1)(g, 2),
			       jac);

    curlPhiI = Mapper::curl(eCurlFun(dofI, g * 3),
			    eCurlFun(dofI, g * 3 + 1),
			    eCurlFun(dofI, g * 3 + 2),
			    jac, 1 / det);

    curlPhiJ = Mapper::curl(eCurlFun(dofJ, g * 3),
			    eCurlFun(dofJ, g * 3 + 1),
			    eCurlFun(dofJ, g * 3 + 2),
			    jac, 1 / det);

    integral1 +=
      ((curlPhiI * curlPhiJ) / mu) * fabs(det) * (*gW1)(g);
  }


  // Loop over Integration Point (Term 2) //
  for(int g = 0; g < G2; g++){
    det = celement.getJacobian((*gC2)(g, 0),
			       (*gC2)(g, 1),
			       (*gC2)(g, 2),
			       invJac);

    invJac.invertInPlace();

    phiI = Mapper::grad(eFun(dofI, g * 3),
			eFun(dofI, g * 3 + 1),
			eFun(dofI, g * 3 + 2),
			invJac);

    phiJ = Mapper::grad(eFun(dofJ, g * 3),
			eFun(dofJ, g * 3 + 1),
			eFun(dofJ, g * 3 + 2),
			invJac);

    integral2 +=
      ((phiI * phiJ) * eps * kSquare) * fabs(det) * (*gW2)(g);
  }

  return integral1 - integral2;
}
