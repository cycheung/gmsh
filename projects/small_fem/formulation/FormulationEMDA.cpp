#include "GroupOfJacobian.h"
#include "GroupOfElement.h"
#include "Quadrature.h"

#include "FormulationEMDA.h"

using namespace std;

FormulationEMDA::
FormulationEMDA(const FunctionSpaceScalar& fs,
                const std::map<Dof, std::complex<double> >& ddmDof){
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

  // DDM //
  this->ddmDof = &ddmDof;
}

FormulationEMDA::~FormulationEMDA(void){
  delete gC;
  delete gW;
  delete jac;
}

complex<double> FormulationEMDA::weak(size_t dofI, size_t dofJ,
                                      size_t elementId) const{
  return complex<double>(0, 0);
}

complex<double> FormulationEMDA::rhs(size_t equationI, size_t elementId) const{
  // Init //
  double phi;
  double det;

  double pxyz[3];
  fullVector<double> xyz(3);

  complex<double> ddmValue;
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

    // Compute ddmValue in the *physical* coordinate
    basis->getReferenceSpace().mapFromABCtoXYZ(element,
                                               (*gC)(g, 0),
                                               (*gC)(g, 1),
                                               (*gC)(g, 2),
                                               pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    ddmValue = interpolate(element, xyz);

    // Integrate
    integral += ddmValue * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}

std::complex<double> FormulationEMDA::
interpolate(const MElement& element, const fullVector<double>& xyz) const{
  // Get Dofs associated to element //
  const vector<Dof>  dof = fspace->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, complex<double> >::const_iterator end = ddmDof->end();
  map<Dof, complex<double> >::const_iterator it;
  vector<double> realCoef(nDof);
  vector<double> imagCoef(nDof);

  for(size_t i = 0; i < nDof; i++){
    it = ddmDof->find(dof[i]);
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

bool FormulationEMDA::isGeneral(void) const{
  return false;
}

complex<double>  FormulationEMDA::weakB(size_t dofI, size_t dofJ,
                                        size_t elementId) const{
  return complex<double>(0, 0);
}

const FunctionSpace& FormulationEMDA::fs(void) const{
  return *fspace;
}
