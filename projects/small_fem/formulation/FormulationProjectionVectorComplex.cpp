#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "GroupOfElement.h"
#include "Quadrature.h"
#include "Mapper.h"

#include "FormulationProjectionVector.h"

using namespace std;

template<>
FormulationProjectionVector<complex<double> >::
FormulationProjectionVector(fullVector<complex<double> >
                                      (*f)(fullVector<double>& xyz),
                            FunctionSpaceVector& fs){

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
  jac = new GroupOfJacobian(*goe, *basis, *gC, "invert");

  // f //
  this->f = f;
}

template<>
FormulationProjectionVector<complex<double> >::
~FormulationProjectionVector(void){
  delete gC;
  delete gW;
  delete jac;
}

template<>
complex<double> FormulationProjectionVector<complex<double> >::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  // Init //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);

  double det;
  double integral = 0;

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians /
  const fullMatrix<double>* myJac;

  const vector<const pair<const fullMatrix<double>*, double>*>& allJac =
    jac->getJacobian(elementId).getInvertJacobianMatrix();

  // Loop over Integration Point //
  const size_t G = gW->size();

  for(size_t g = 0; g < G; g++){
    myJac = allJac[g]->first;
    det   = allJac[g]->second;

    Mapper::hCurl(eFun, dofI, g, *myJac, phiI);
    Mapper::hCurl(eFun, dofJ, g, *myJac, phiJ);

    integral += (phiI * phiJ) * fabs(det) * (*gW)(g);
  }

  return complex<double>(integral, 0);
}

template<>
complex<double> FormulationProjectionVector<complex<double> >::
rhs(size_t equationI, size_t elementId) const{

  // Init //
  fullVector<double> phi(3);

  fullVector<double>           xyz(3);
  double                       pxyz[3];
  fullVector<complex<double> > fxyz;

  double          det;
  complex<double> integral = complex<double>(0, 0);
  complex<double> tmp;

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians /
  const fullMatrix<double>* myJac;

  const vector<const pair<const fullMatrix<double>*, double>*>& allJac =
    jac->getJacobian(elementId).getInvertJacobianMatrix();

  // Loop over Integration Point //
  const size_t G = gW->size();

  for(size_t g = 0; g < G; g++){
    // Jacobian
    myJac = allJac[g]->first;
    det   = allJac[g]->second;

    // Basis Function
    Mapper::hCurl(eFun, equationI, g, *myJac, phi);

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
    tmp = fxyz(0) * phi(0) + fxyz(1) * phi(1) + fxyz(2) * phi(2);
    integral += tmp * fabs(det) * (*gW)(g);
  }

  return integral;
}

template<>
bool FormulationProjectionVector<complex<double> >::isGeneral(void) const{

  return false;
}

template<>
complex<double> FormulationProjectionVector<complex<double> >::
weakB(size_t dofI, size_t dofJ, size_t elementId) const{

  return complex<double>(0, 0);
}

template<>
const FunctionSpace& FormulationProjectionVector<complex<double> >::
fs(void) const{

  return *fspace;
}
