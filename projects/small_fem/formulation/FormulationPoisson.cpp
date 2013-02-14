#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Jacobian.h"

#include "Exception.h"
#include "FormulationPoisson.h"

using namespace std;

// Source Terms //
const unsigned int FormulationPoisson::sourceOrder = 1;

double FormulationPoisson::gSource(fullVector<double>& xyz){
  return 1;
}

// Poisson //
FormulationPoisson::FormulationPoisson(GroupOfElement& goe,
				       unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Poisson formulation of order 0");

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature Data (LHS) //
  // NB: We need to integrad a grad * grad !
  //     and order(grad f) = order(f) - 1
  fullMatrix<double> gCL;
  fullVector<double> gWL;

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), 2 * (order - 1), gCL, gWL);

  // Gaussian Quadrature Data (RHS) //
  // NB: We need to integrad a f * gSource !
  //     and order(gSource) = sourceOrder
  fullMatrix<double> gCR;
  fullVector<double> gWR;

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), order + sourceOrder, gCR, gWR);

  // Local Terms //
  basis->preEvaluateDerivatives(gCL);
  basis->preEvaluateFunctions(gCR);
  goe.orientAllElements(*basis);

  Jacobian jacL(goe, gCL);
  Jacobian jacR(goe, gCR);
  jacL.computeInvertJacobians();
  jacR.computeJacobians();

  localTermsL = new TermGradGrad(jacL, *basis, gWL);
  localTermsR = new TermProjectionField(jacR, *basis, gWR, gCR, gSource);
}

FormulationPoisson::~FormulationPoisson(void){
  delete basis;
  delete fspace;

  delete localTermsL;
  delete localTermsR;
}

double FormulationPoisson::weak(unsigned int dofI, unsigned int dofJ,
                                unsigned int elementId) const{

  return localTermsL->getTerm(dofI, dofJ, elementId);
}

double FormulationPoisson::rhs(unsigned int equationI,
                               unsigned int elementId) const{

  return localTermsR->getTerm(0, equationI, elementId);
}

bool FormulationPoisson::isGeneral(void) const{
  return false;
}

double FormulationPoisson::weakB(unsigned int dofI,
                                 unsigned int dofJ,
                                 unsigned int elementId) const{
  return 0;
}

const FunctionSpace& FormulationPoisson::fs(void) const{
  return *fspace;
}
