#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Mapper.h"

#include "Exception.h"
#include "FormulationUpdateOO2.h"

using namespace std;

FormulationUpdateOO2::
FormulationUpdateOO2(const FunctionSpaceScalar& fs,
                     std::complex<double> a,
                     std::complex<double> b,
                     const std::map<Dof, std::complex<double> >& solution,
                     const std::map<Dof, std::complex<double> >& oldG){
  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(0);

  // a & b //
  this->a = a;
  this->b = b;

  // Domain //
  this->goe = &fs.getSupport();

  // Gaussian Quadrature (Field - Field) //
  Quadrature gaussFF(goe->get(0).getType(), basis->getOrder(), 2);

  gCFF = new fullMatrix<double>(gaussFF.getPoints());
  gWFF = new fullVector<double>(gaussFF.getWeights());

  // Gaussian Quadrature (Grad - Grad) //
  Quadrature gaussGG(goe->get(0).getType(), basis->getOrder() - 1, 2);

  gCGG = new fullMatrix<double>(gaussGG.getPoints());
  gWGG = new fullVector<double>(gaussGG.getWeights());

  // Pre-evalution //
  basis->preEvaluateFunctions(*gCFF);
  basis->preEvaluateDerivatives(*gCGG);

  jacFF = new GroupOfJacobian(*goe, *basis, *gCFF, "jacobian");
  jacGG = new GroupOfJacobian(*goe, *basis, *gCGG, "invert");

  // DDM //
  this->solution = &solution;
  this->oldG     = &oldG;
}

FormulationUpdateOO2::~FormulationUpdateOO2(void){
  delete gCFF;
  delete gWFF;
  delete jacFF;

  delete gCGG;
  delete gWGG;
  delete jacGG;
}

complex<double> FormulationUpdateOO2::
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
    jacFF->getJacobian(elementId).getJacobianMatrix();

  // Loop over Integration Point //
  const size_t G = gWFF->size();

  for(size_t g = 0; g < G; g++){
    det   = allJac[g]->second;

    phiI = eFun(dofI, g);
    phiJ = eFun(dofJ, g);

    integral += phiI * phiJ * fabs(det) * (*gWFF)(g);
  }

  return complex<double>(integral, 0);
}

complex<double> FormulationUpdateOO2::
rhs(size_t equationI, size_t elementId) const{
  // Init //
  size_t G;

  double det;
  const fullMatrix<double>* jac;

  double phi;
  fullVector<double> gradPhi(3);

  double pxyz[3];
  fullVector<double> xyz(3);
  complex<double> oldGValue;
  complex<double> solutionValue;
  fullVector<complex<double> > gradValue;

  complex<double> sub;

  complex<double> integral = complex<double>(0, 0);

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eFunFF =
    basis->getPreEvaluatedFunctions(element);

  // Get Grad Basis Functions //
  const fullMatrix<double>& eFunGG =
    basis->getPreEvaluatedDerivatives(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& allJacFF =
    jacFF->getJacobian(elementId).getJacobianMatrix();

  const vector<const pair<const fullMatrix<double>*, double>*>& allJacGG =
    jacGG->getJacobian(elementId).getInvertJacobianMatrix();

  // Loop over Integration Point (Field - Field) //
  G = gWFF->size();

  for(size_t g = 0; g < G; g++){
    // Compute phi
    det = allJacFF[g]->second;
    phi = eFunFF(equationI, g);

    // Get *physical* coordinate
    basis->getReferenceSpace().mapFromABCtoXYZ(element,
                                               (*gCFF)(g, 0),
                                               (*gCFF)(g, 1),
                                               (*gCFF)(g, 2),
                                               pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    // OldG & solution in *physical* coordinate
    oldGValue     = interpolate(element, xyz, *oldG);
    solutionValue = interpolate(element, xyz, *solution);
    sub           = complex<double>(2, 0) * a * solutionValue - oldGValue;

    // Integrate
    integral += sub * phi * fabs(det) * (*gWFF)(g);
  }

  // Loop over Integration Point (Grad - Grad) //
  G = gWGG->size();

  for(size_t g = 0; g < G; g++){
    // Compute gradPhi
    det = allJacGG[g]->second;
    jac = allJacGG[g]->first;

    Mapper::hCurl(eFunGG, equationI, g, *jac, gradPhi);

    // Get *physical* coordinate
    basis->getReferenceSpace().mapFromABCtoXYZ(element,
                                               (*gCGG)(g, 0),
                                               (*gCGG)(g, 1),
                                               (*gCGG)(g, 2),
                                               pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    // grad(solution) in *physical* coordinate
    gradValue = interpolateGrad(element, xyz, *solution);

    // Integrate
    integral +=
      complex<double>(-2, 0) * b *
      (gradValue(0) * complex<double>(gradPhi(0), 0) +
       gradValue(1) * complex<double>(gradPhi(1), 0) +
       gradValue(2) * complex<double>(gradPhi(2), 0)) * fabs(det) * (*gWGG)(g);
  }

  return integral;
}

std::complex<double> FormulationUpdateOO2::
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

fullVector<std::complex<double> > FormulationUpdateOO2::
interpolateGrad(const MElement& element,
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
  fullVector<double> re = fspace->interpolateDerivative(element, realCoef, xyz);
  fullVector<double> im = fspace->interpolateDerivative(element, imagCoef, xyz);

  // Return //
  if(re.size() != 3)
    throw Exception("Snif");

  fullVector<complex<double> > ret(3);

  ret(0) = complex<double>(re(0), im(0));
  ret(1) = complex<double>(re(1), im(1));
  ret(2) = complex<double>(re(2), im(2));

  return ret;
}

bool FormulationUpdateOO2::isGeneral(void) const{
  return false;
}

complex<double>  FormulationUpdateOO2::weakB(size_t dofI, size_t dofJ,
                                             size_t elementId) const{
  return complex<double>(0, 0);
}

const FunctionSpace& FormulationUpdateOO2::fs(void) const{
  return *fspace;
}
