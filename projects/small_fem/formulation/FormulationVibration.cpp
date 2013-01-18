#include <cmath>

#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Mapper.h"

#include "Exception.h"
#include "FormulationVibration.h"

using namespace std;

FormulationVibration::FormulationVibration(const GroupOfElement& goe,
					   unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Vibration formulation of order 0");

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature Data (LHS) //
  // NB: We need to integrad a grad * grad !
  //     and order(grad f) = order(f) - 1
  gCL = new fullMatrix<double>();
  gWL = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), (order - 1) + (order - 1), *gCL, *gWL);

  GL = gWL->size(); // Nbr of Gauss points

  // PreEvaluate
  basis->preEvaluateDerivatives(*gCL);
}

FormulationVibration::~FormulationVibration(void){
  delete gCL;
  delete gWL;
  delete basis;
  delete fspace;
}

double FormulationVibration::weakA(int dofI, int dofJ,
				   const GroupOfDof& god) const{

  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double> invJac(3, 3);
  double integral = 0;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedDerivatives(element);

  // Loop over Integration Point //
  for(int g = 0; g < GL; g++){
    double det = celement.getJacobian((*gCL)(g, 0),
				      (*gCL)(g, 1),
				      (*gCL)(g, 2),
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

    integral += phiI * phiJ * fabs(det) * (*gWL)(g);
  }

  return integral;
}
