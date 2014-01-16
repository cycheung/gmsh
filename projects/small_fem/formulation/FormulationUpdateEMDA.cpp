#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationUpdateEMDA.h"

using namespace std;

FormulationUpdateEMDA::
FormulationUpdateEMDA(const FunctionSpaceScalar& fs,
                      double k,
                      double chi,
                      const std::map<Dof, std::complex<double> >& solution,
                      const std::map<Dof, std::complex<double> >& oldG){
  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(0);

  // Wavenumber & Chi //
  this->k   = k;
  this->chi = chi;

  // Domain //
  this->goe = &fs.getSupport();

  // Gaussian Quadrature //
  Quadrature gauss(goe->get(0).getType(), basis->getOrder(), 2);

  gC = new fullMatrix<double>(gauss.getPoints());
  gW = new fullVector<double>(gauss.getWeights());

  // Pre-evalution //
  basis->preEvaluateFunctions(*gC);
  jac = new GroupOfJacobian(*goe, *basis, *gC, "jacobian");

  // DDM //
  this->solution = &solution;
  this->oldG     = &oldG;
}

FormulationUpdateEMDA::~FormulationUpdateEMDA(void){
  delete gC;
  delete gW;
  delete jac;
}

complex<double> FormulationUpdateEMDA::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
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

complex<double> FormulationUpdateEMDA::
rhs(size_t equationI, size_t elementId) const{
  // Init //
  double phi;
  double det;

  double pxyz[3];
  fullVector<double> xyz(3);
  complex<double> oldGValue;
  complex<double> solutionValue;
  complex<double> sub;

  complex<double> integral = complex<double>(0, 0);

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

    // Get *physical* coordinate
    basis->getReferenceSpace().mapFromABCtoXYZ(element,
                                               (*gC)(g, 0),
                                               (*gC)(g, 1),
                                               (*gC)(g, 2),
                                               pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    // OldG & solution in *physical* coordinate
    oldGValue     = interpolate(element, xyz, *oldG);
    solutionValue = interpolate(element, xyz, *solution);
    sub           =
      complex<double>(2 * chi, -2 * k) * solutionValue - oldGValue;

    // Integrate
    integral += sub * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}

std::complex<double> FormulationUpdateEMDA::
interpolate(const MElement& element,
            const fullVector<double>& xyz,
            const std::map<Dof, std::complex<double> >& f) const{

  // Get Dofs associated to element //
  const vector<Dof>  dof = fspace->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, complex<double> >::const_iterator end = f.end();
  map<Dof, complex<double> >::const_iterator it;
  vector<double> realCoef(nDof);
  vector<double> imagCoef(nDof);

  for(size_t i = 0; i < nDof; i++){
    it = f.find(dof[i]);
    if(it == end)
      throw Exception("Snif");

    realCoef[i] = it->second.real();
    imagCoef[i] = it->second.imag();
  }

  // Interpolate
  double real = fspace->interpolate(element, realCoef, xyz);
  double imag = fspace->interpolate(element, imagCoef, xyz);

  return complex<double>(real, imag);
}

bool FormulationUpdateEMDA::isGeneral(void) const{
  return false;
}

complex<double>  FormulationUpdateEMDA::weakB(size_t dofI, size_t dofJ,
                                              size_t elementId) const{
  return complex<double>(0, 0);
}

const FunctionSpace& FormulationUpdateEMDA::fs(void) const{
  return *fspace;
}
