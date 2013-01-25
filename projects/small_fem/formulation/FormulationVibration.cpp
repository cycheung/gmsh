#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Jacobian.h"

#include "Exception.h"
#include "FormulationVibration.h"

using namespace std;

FormulationVibration::FormulationVibration(GroupOfElement& goe,
					   unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Vibration formulation of order 0");

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature Data //
  // NB: We need to integrad a grad * grad !
  //     and order(grad f) = order(f) - 1
  fullMatrix<double> gC;
  fullVector<double> gW;

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), (order - 1) + (order - 1), gC, gW);

  // Local Terms //
  basis->preEvaluateDerivatives(gC);
  goe.orientAllElements(*basis);

  Jacobian jac(goe, gC);
  jac.computeInvertJacobians();

  localTerms = new TermHCurl(jac, *basis, gW);
}

FormulationVibration::~FormulationVibration(void){
  delete basis;
  delete fspace;

  delete localTerms;
}

double FormulationVibration::weakA(unsigned int dofI, unsigned int dofJ,
				   const GroupOfDof& god) const{

  return localTerms->getTerm(dofI, dofJ, god);
}

double FormulationVibration::weakB(unsigned int dofI, unsigned int dofJ,
					  const GroupOfDof& god) const{
  throw
    Exception
    ("Vibration is a Non General Eigenvalue problem, and don't need a B matrix");
}

bool FormulationVibration::isGeneral(void) const{
  return false;
}

const FunctionSpace& FormulationVibration::fs(void) const{
  return *fspace;
}
