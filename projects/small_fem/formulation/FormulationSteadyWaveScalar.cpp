#include "Exception.h"
#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "Polynomial.h"
#include "Mapper.h"

#include "FormulationSteadyWaveScalar.h"

using namespace std;

// Pi  = atan(1) * 4
// Mu  = 4 * Pi * 10^-7
// Eps = 8.85418781762 * 10^−12
//const double FormulationSteadyWaveScalar::mu  = 4 * atan(1) * 4 * 1E-7;
//const double FormulationSteadyWaveScalar::eps = 8.85418781762E-12;

const double FormulationSteadyWaveScalar::mu  = 1;
const double FormulationSteadyWaveScalar::eps = 1;

FormulationSteadyWaveScalar::FormulationSteadyWaveScalar(const GroupOfElement& goe,
							 double k,
							 unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Scalar SteadyWave formulation of order 0");

  // Wave Number Squared //
  kSquare = k * k;

  // Function Space & Basis//
  fspace = new FunctionSpaceNode(goe, order);
  basis  = &fspace->getBasis(0);

  // Gaussian Quadrature Data (Term One) //
  // NB: We need to integrad a grad * grad !
  //     and order(rot f) = order(f) - 1
  gC1 = new fullMatrix<double>();
  gW1 = new fullVector<double>();

  // Gaussian Quadrature Data (Term Two) //
  // NB: We need to integrad a f * f !
  gC2 = new fullMatrix<double>();
  gW2 = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), 2 * (order - 1), *gC1, *gW1);
  gaussIntegration::get(goe.get(0).getType(), 2 *  order     , *gC2, *gW2);

  // Nbr of Gauss points
  G1 = gW1->size();
  G2 = gW2->size();

  // PreEvaluate
  basis->preEvaluateDerivatives(*gC1);
  basis->preEvaluateFunctions(*gC2);
}

FormulationSteadyWaveScalar::~FormulationSteadyWaveScalar(void){
  delete gC1;
  delete gW1;
  delete gC2;
  delete gW2;
  delete fspace;
}

double FormulationSteadyWaveScalar::weak(int dofI, int dofJ,
					 const GroupOfDof& god) const{
  // Init Some Stuff //
  fullVector<double> gradPhiI(3);
  fullVector<double> gradPhiJ(3);
  double phiI;
  double phiJ;

  fullMatrix<double> invJac(3, 3);

  double integral1 = 0;
  double integral2 = 0;
  double det;

  // Get Element and Basis Functions (+ Grad) //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);

  const fullMatrix<double>& eGradFun =
    basis->getPreEvaluatedDerivatives(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Loop over Integration Point (Term 1) //
  for(int g = 0; g < G1; g++){
    det = celement.getJacobian((*gC1)(g, 0),
			       (*gC1)(g, 1),
			       (*gC1)(g, 2),
			       invJac);

    invJac.invertInPlace();

    gradPhiI = Mapper::grad(eGradFun(dofI, g * 3),
			    eGradFun(dofI, g * 3 + 1),
			    eGradFun(dofI, g * 3 + 2),
			    invJac);

    gradPhiJ = Mapper::grad(eGradFun(dofJ, g * 3),
			    eGradFun(dofJ, g * 3 + 1),
			    eGradFun(dofJ, g * 3 + 2),
			    invJac);

    integral1 +=
      ((gradPhiI * gradPhiJ) / mu) * fabs(det) * (*gW1)(g);
  }


  // Loop over Integration Point (Term 2) //
  for(int g = 0; g < G2; g++){
    det = celement.getJacobian((*gC2)(g, 0),
			       (*gC2)(g, 1),
			       (*gC2)(g, 2),
			       invJac);


    phiI = eFun(dofI, g);
    phiJ = eFun(dofJ, g);

    integral2 +=
      ((phiI * phiJ) * eps * kSquare) * fabs(det) * (*gW2)(g);
  }

  return integral1 - integral2;
}
