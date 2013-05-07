#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationEigenFrequency.h"

using namespace std;

// Pi  = atan(1) * 4
// Mu  = 4 * Pi * 10^-7
// Eps = 8.85418781762 * 10^−12
//const double FormulationEigenFrequency::mu  = 4 * atan(1) * 4 * 1E-7;
//const double FormulationEigenFrequency::eps = 8.85418781762E-12;

const double FormulationEigenFrequency::mu  = 1;
const double FormulationEigenFrequency::eps = 1;

FormulationEigenFrequency::FormulationEigenFrequency(GroupOfElement& goe,
						     unsigned int order){
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

  GroupOfJacobian jac1(goe, gC1, "jacobian");
  GroupOfJacobian jac2(goe, gC2, "invert");

  localTerms1 = new TermCurlCurl(jac1, *basis, gW1);
  localTerms2 = new TermGradGrad(jac2, *basis, gW2);
}

FormulationEigenFrequency::~FormulationEigenFrequency(void){
  delete basis;
  delete fspace;

  delete localTerms1;
  delete localTerms2;
}

double FormulationEigenFrequency::weak(unsigned int dofI, unsigned int dofJ,
                                       unsigned int elementId) const{

  return localTerms1->getTerm(dofI, dofJ, elementId) / mu;
}

double FormulationEigenFrequency::weakB(unsigned int dofI, unsigned int dofJ,
					unsigned int elementId) const{

  return localTerms2->getTerm(dofI, dofJ, elementId) * eps;
}

double FormulationEigenFrequency::rhs(unsigned int dofI,
                                      unsigned int elementId) const{
  return 0;
}

bool FormulationEigenFrequency::isGeneral(void) const{
  return true;
}

const FunctionSpace& FormulationEigenFrequency::fs(void) const{
  return *fspace;
}
