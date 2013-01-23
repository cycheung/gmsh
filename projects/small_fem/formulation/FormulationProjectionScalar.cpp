#include <cmath>

#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Mapper.h"

#include "FormulationProjectionScalar.h"

using namespace std;

FormulationProjectionScalar::
FormulationProjectionScalar(double (*f)(fullVector<double>& xyz),
			    FunctionSpaceScalar& fs){
  // Save f //
  this->f = f;

  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(0);
  goe    = &fs.getSupport();

  // Gaussian Quadrature Data  //
  // NB: We need to integrad f_i * f_j or f_i * g
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe->get(0).getType(), 2 * basis->getOrder(), *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // PreEvaluate
  basis->preEvaluateFunctions(*gC);

  // Fast Assembly //
  goe->orientAllElements(*basis);
  jac = new Jacobian(*goe, *gC);
  jac->computeJacobians();

  orientationStat = &goe->getOrientationStats();
  nOrientation    = basis->getNOrientation();
  nFunction       = basis->getNFunction();
  nElement        = goe->getNumber();

  eMap = new map<const MElement*, pair<unsigned int, unsigned int> >;

  computeC();
  computeB();
  computeA();

  deleteCB();
  delete jac;
}

FormulationProjectionScalar::~FormulationProjectionScalar(void){
  delete gC;
  delete gW;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete laM[s];

  delete[] laM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete raM[s];

  delete[] raM;

  delete eMap;
}

double FormulationProjectionScalar::weak(int dofI, int dofJ,
                                         const GroupOfDof& god) const{

  map<const MElement*, pair<unsigned int, unsigned int> >::iterator
    index = eMap->find(&god.getGeoElement());

  return (*laM[index->second.first])
    (index->second.second, dofI * nFunction + dofJ);
}

double FormulationProjectionScalar::rhs(int equationI,
					const GroupOfDof& god) const{

  map<const MElement*, pair<unsigned int, unsigned int> >::iterator
    index = eMap->find(&god.getGeoElement());

  return (*raM[index->second.first])(index->second.second, equationI);
}

void FormulationProjectionScalar::computeC(void){
  unsigned int l;

  // Alloc //
  lcM = new fullMatrix<double>*[nOrientation];
  rcM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    lcM[s] = new fullMatrix<double>(G, nFunction * nFunction);

  for(unsigned int s = 0; s < nOrientation; s++)
    rcM[s] = new fullMatrix<double>(G, nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis->getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(int g = 0; g < G; g++){
      // Reset Function Counter
      l = 0;

      // Loop on Functions
      for(unsigned int i = 0; i < nFunction; i++)
        (*rcM[s])(g, i) = (*gW)(g) * phi(i, g);

      for(unsigned int i = 0; i < nFunction; i++){
        for(unsigned int j = 0; j < nFunction; j++){
          (*lcM[s])(g, l) = (*rcM[s])(g, i) * phi(j, g);
          l++;
        }
      }
    }
  }
}

void FormulationProjectionScalar::computeB(void){
  const vector<const MElement*>& element = goe->getAll();
  unsigned int                   offset  = 0;
  unsigned int j;

  fullVector<double> xyz(3);
  SPoint3            pxyz;
  double             fxyz;

  // Alloc //
  lbM = new fullMatrix<double>*[nOrientation];
  rbM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    lbM[s] = new fullMatrix<double>((*orientationStat)[s], G);

  for(unsigned int s = 0; s < nOrientation; s++)
    rbM[s] = new fullMatrix<double>((*orientationStat)[s], G);


  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Reset Element Counter
    j = 0;

    // Loop on Elements
    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      // Add to eMap
      eMap->insert(pair<const MElement*, pair<unsigned int, unsigned int> >
                     (element[e], pair<unsigned int, unsigned int>(s, j)));

      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& jacM =
        jac->getJacobian(*element[e]);

      // Loop on Gauss Points
      for(int g = 0; g < G; g++)
        (*lbM[s])(j, g) = fabs(jacM[g]->second);

      for(int g = 0; g < G; g++){
        // Compute f in the *physical* coordinate
        const_cast<MElement*>(element[e])
          ->pnt((*gC)(g, 0),
                (*gC)(g, 1),
                (*gC)(g, 2),
                pxyz);

        xyz(0) = pxyz.x();
        xyz(1) = pxyz.y();
        xyz(2) = pxyz.z();

        fxyz = f(xyz);

        (*rbM[s])(j, g) = (*lbM[s])(j, g) * fxyz;
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}

void FormulationProjectionScalar::computeA(void){
  // Alloc //
  laM = new fullMatrix<double>*[nOrientation];
  raM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    laM[s] = new fullMatrix<double>((*orientationStat)[s], nFunction * nFunction);

  for(unsigned int s = 0; s < nOrientation; s++)
    raM[s] = new fullMatrix<double>((*orientationStat)[s], nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // GEMM doesn't like matrices with 0 Elements
    if((*orientationStat)[s]){
      laM[s]->gemm(*lbM[s], *lcM[s]);
      raM[s]->gemm(*rbM[s], *rcM[s]);
    }
  }
}

void FormulationProjectionScalar::deleteCB(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete lcM[s];

  delete[] lcM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete lbM[s];

  delete[] lbM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete rcM[s];

  delete[] rcM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete rbM[s];

  delete[] rbM;
}

/*
double FormulationProjectionScalar::weak(int dofI, int dofJ,
				   const GroupOfDof& god) const{
  // Init //
  double det;
  double phiI;
  double phiJ;
  // fullMatrix<double> jac(3, 3);
  double integral = 0;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  // MElement&      celement = const_cast<MElement&>(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  const vector<const pair<const fullMatrix<double>*, double>*>& jacM =
    jac->getJacobian(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    det = jacM[g]->second;

    // det = celement.getJacobian((*gC)(g, 0),
    //                            (*gC)(g, 1),
    //			          (*gC)(g, 2),
    //			          jac);

    phiI = eFun(dofI, g);
    phiJ = eFun(dofJ, g);

    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return integral;
}
*/
/*
double FormulationProjectionScalar::rhs(int equationI,
					const GroupOfDof& god) const{
  // Init //
  double phi;
  double det;

  fullVector<double> xyz(3);
  SPoint3            pxyz;
  double             fxyz;

  double integral = 0;
  // fullMatrix<double> jac(3, 3);

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  const vector<const pair<const fullMatrix<double>*, double>*>& jacM =
    jac->getJacobian(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    // Compute phi
    det = jacM[g]->second;

    // det = celement.getJacobian((*gC)(g, 0),
    //                            (*gC)(g, 1),
    //                            (*gC)(g, 2),
    //                            jac);

    phi = eFun(equationI, g);

    // Compute f in the *physical* coordinate
    celement.pnt((*gC)(g, 0),
		 (*gC)(g, 1),
		 (*gC)(g, 2),
		 pxyz);

    xyz(0) = pxyz.x();
    xyz(1) = pxyz.y();
    xyz(2) = pxyz.z();

    fxyz = f(xyz);

    // Integrate
    integral += fxyz * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}
*/
