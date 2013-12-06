#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationNeumann.h"

using namespace std;

FormulationNeumann::
FormulationNeumann(GroupOfElement& goe, double omega, size_t order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Scalar SteadyWave formulation of order 0");

  // Pulsation Squared //
  this->omega = omega;

  // Function Space & Basis//
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature //
  Quadrature gaussFF(goe.get(0).getType(), order, 2);
  const fullMatrix<double>& gC = gaussFF.getPoints();
  const fullVector<double>& gW = gaussFF.getWeights();

  // Local Terms //
  basis->preEvaluateFunctions(gC);

  GroupOfJacobian jac(goe, *basis, gC, "jacobian");

  localTerms = new TermFieldField(jac, *basis, gW);
}

FormulationNeumann::~FormulationNeumann(void){
  delete basis;
  delete fspace;

  delete localTerms;
}

complex<double> FormulationNeumann::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{
  return
    complex<double>(0, -1 * omega * localTerms->getTerm(dofI, dofJ, elementId));
}

complex<double> FormulationNeumann::rhs(size_t equationI,
                                        size_t elementId) const{
  return complex<double>(0, 0);
}

bool FormulationNeumann::isGeneral(void) const{
  return false;
}

complex<double> FormulationNeumann::weakB(size_t dofI, size_t dofJ,
                                          size_t elementId) const{
  return complex<double>(0, 0);
}


const FunctionSpace& FormulationNeumann::fs(void) const{
  return *fspace;
}
