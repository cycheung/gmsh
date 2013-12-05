#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "GroupOfElement.h"
#include "Quadrature.h"

#include "FormulationProjectionScalar.h"

using namespace std;

template<>
FormulationProjectionScalar<complex<double> >::
FormulationProjectionScalar(complex<double> (*f)(fullVector<double>& xyz),
                            FunctionSpaceScalar& fs){
  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(0);

  // Domain //
  goe = &fs.getSupport();

  // Gaussian Quadrature //
  Quadrature gauss(goe->get(0).getType(), basis->getOrder(), 2);

  gC = new fullMatrix<double>(gauss.getPoints());
  gW = new fullVector<double>(gauss.getWeights());

  // Pre-evalution //
  basis->preEvaluateFunctions(*gC);
  jac = new GroupOfJacobian(*goe, *basis, *gC, "jacobian");

  // f //
  this->f = f;
}

template<>
FormulationProjectionScalar<complex<double> >::
~FormulationProjectionScalar(void){
  delete gC;
  delete gW;
  delete jac;
}

template<>
complex<double> FormulationProjectionScalar<complex<double> >::
weak(size_t dofI, size_t dofJ,size_t elementId) const{

  // Init //
  double phiI;
  double phiJ;

  double det;
  double integral = 0;

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& allJac =
    jac->getJacobian(elementId).getJacobianMatrix();

  // Loop over Integration Point //
  const size_t G = gW->size();

  for(size_t g = 0; g < G; g++){
    det   = allJac[g]->second;

    phiI = eFun(dofI, g);
    phiJ = eFun(dofJ, g);

    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return complex<double>(integral, 0);
}

template<>
complex<double> FormulationProjectionScalar<complex<double> >::
rhs(size_t equationI, size_t elementId) const{

  // Init //
  double phi;
  double det;

  fullVector<double> xyz(3);
  double             pxyz[3];
  complex<double>    fxyz;

  complex<double>    integral = complex<double>(0, 0);

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& allJac =
    jac->getJacobian(elementId).getJacobianMatrix();

  // Loop over Integration Point //
  const size_t G = gW->size();

  for(size_t g = 0; g < G; g++){
    // Compute phi
    det = allJac[g]->second;
    phi = eFun(equationI, g);

    // Compute f in the *physical* coordinate
    basis->getReferenceSpace().mapFromABCtoXYZ(element,
                                               (*gC)(g, 0),
                                               (*gC)(g, 1),
                                               (*gC)(g, 2),
                                               pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    fxyz = f(xyz);

    // Integrate
    integral += fxyz * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}

template<>
bool FormulationProjectionScalar<complex<double> >::isGeneral(void) const{

  return false;
}

template<>
complex<double>  FormulationProjectionScalar<complex<double> >::
weakB(size_t dofI, size_t dofJ, size_t elementId) const{

  return complex<double>(0, 0);
}

template<>
const FunctionSpace& FormulationProjectionScalar<complex<double> >::
fs(void) const{

  return *fspace;
}
