#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Jacobian.h"

#include "Exception.h"
#include "FormulationLaplace.h"

using namespace std;

FormulationLaplace::FormulationLaplace(GroupOfElement& goe,
				       unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Laplace formulation of order 0");

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature Data //
  fullMatrix<double> gC;
  fullVector<double> gW;

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), 2 * (order - 1), gC, gW);

  // Local Terms //
  basis->preEvaluateDerivatives(gC);
  goe.orientAllElements(*basis);

  Jacobian jac(goe, gC);
  jac.computeInvertJacobians();

  localTerms = new TermGradGrad(jac, *basis, gW);
}

FormulationLaplace::~FormulationLaplace(void){
  delete basis;
  delete fspace;
  delete localTerms;
}

double FormulationLaplace::weak(unsigned int dofI, unsigned int dofJ,
				const GroupOfDof& god) const{

  return localTerms->getTerm(dofI, dofJ, god);
}


double FormulationLaplace::rhs(unsigned int equationI,
                               const GroupOfDof& god) const{
  return 0;
}

bool FormulationLaplace::isGeneral(void) const{
  return false;
}

double FormulationLaplace::weakB(unsigned int dofI,
                                 unsigned int dofJ,
                                 const GroupOfDof& god) const{
  return 0;
}

const FunctionSpace& FormulationLaplace::fs(void) const{
  return *fspace;
}
