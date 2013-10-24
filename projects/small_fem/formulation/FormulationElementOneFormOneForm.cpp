#include "FormulationElementOneFormOneForm.h"

#include "BasisGenerator.h"
#include "Quadrature.h"
#include "Exception.h"

FormulationElementOneFormOneForm::
FormulationElementOneFormOneForm(GroupOfElement& goe,
                                 size_t order,
                                 size_t orientation){
  // Basis //
  basis = BasisGenerator::generate(goe.get(0).getType(),
                                   1, order, "hierarchical");

  fspace = new FunctionSpaceVector(goe, *basis); // Here for compatibility

  // Save orientation //
  this->orientation = orientation;

  // Gaussian Quadrature //
  Quadrature gauss(goe.get(0).getType(), order - 1, 2);

  gC = new fullMatrix<double>(gauss.getPoints());
  gW = new fullVector<double>(gauss.getWeights());

  // Local Terms //
  basis->preEvaluateFunctions(*gC);
}

FormulationElementOneFormOneForm::~FormulationElementOneFormOneForm(void){
  delete basis;
  delete fspace;
  delete gC;
  delete gW;
}

double FormulationElementOneFormOneForm::weak(size_t dofI, size_t dofJ,
                                              size_t elementId) const{
  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);

  double integral = 0;

  // Get Basis Functions //
  const fullMatrix<double>& fun =
    basis->getPreEvaluatedFunctions(orientation);

  // Loop over Integration Point //
  const int G = gW->size();

  for(int g = 0; g < G; g++){
    phiI(0) = fun(dofI, g * 3);
    phiI(1) = fun(dofI, g * 3 + 1);
    phiI(2) = fun(dofI, g * 3 + 2);

    phiJ(0) = fun(dofJ, g * 3);
    phiJ(1) = fun(dofJ, g * 3 + 1);
    phiJ(2) = fun(dofJ, g * 3 + 2);

    integral += (phiI * phiJ) * (*gW)(g);
  }

  return integral;
}


double FormulationElementOneFormOneForm::rhs(size_t equationI,
                                             size_t elementId) const{
  return 0;
}

bool FormulationElementOneFormOneForm::isGeneral(void) const{
  return false;
}

double FormulationElementOneFormOneForm::weakB(size_t dofI,
                                               size_t dofJ,
                                               size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationElementOneFormOneForm::fs(void) const{
  return *fspace;
}
