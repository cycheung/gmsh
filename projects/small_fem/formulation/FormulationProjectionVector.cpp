#include <cmath>

#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Mapper.h"

#include "FormulationProjectionVector.h"

using namespace std;

FormulationProjectionVector::
FormulationProjectionVector(fullVector<double> (*f)(fullVector<double>& xyz),
			    FunctionSpaceVector& fs){
  // Vector to Project //
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
  jac->computeInvertJacobians();

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

FormulationProjectionVector::~FormulationProjectionVector(void){
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

double FormulationProjectionVector::weak(int dofI, int dofJ,
                                         const GroupOfDof& god) const{

  map<const MElement*, pair<unsigned int, unsigned int> >::iterator
    index = eMap->find(&god.getGeoElement());

  return (*laM[index->second.first])
    (index->second.second, dofI * nFunction + dofJ);
}

double FormulationProjectionVector::rhs(int equationI,
					const GroupOfDof& god) const{

  map<const MElement*, pair<unsigned int, unsigned int> >::iterator
    index = eMap->find(&god.getGeoElement());

  return (*raM[index->second.first])(index->second.second, equationI);
}

void FormulationProjectionVector::computeC(void){
  unsigned int k;
  unsigned int l;

  // Alloc //
  lcM = new fullMatrix<double>*[nOrientation];
  rcM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    lcM[s] = new fullMatrix<double>(9 * G, nFunction * nFunction);

  for(unsigned int s = 0; s < nOrientation; s++)
    rcM[s] = new fullMatrix<double>(3 * G, nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis->getPreEvaluatedFunctions(s);

    // Reset Gauss Point Counter
    k = 0;

    // Loop on Gauss Points
    for(int g = 0; g < G; g++){
      for(unsigned int a = 0; a < 3; a++){
        for(unsigned int b = 0; b < 3; b++){
          // Reset Function Counter
          l = 0;

          // Loop on Functions
          for(unsigned int j = 0; j < nFunction; j++)
            (*rcM[s])(g * 3 + b, j) = (*gW)(g) * phi(j, g * 3 + b);

          for(unsigned int i = 0; i < nFunction; i++){
            for(unsigned int j = 0; j < nFunction; j++){
              (*lcM[s])(k, l) = (*rcM[s])(g * 3 + b, j) * phi(i, g * 3 + a);
              l++;
            }
          }

          k++;
        }
      }
    }
  }
}

void FormulationProjectionVector::computeB(void){
  const vector<const MElement*>& element = goe->getAll();
  unsigned int                   offset  = 0;

  unsigned int j;
  unsigned int k;

  fullVector<double> xyz(3);
  SPoint3            pxyz;
  fullVector<double> fxyz;

  // Alloc //
  lbM = new fullMatrix<double>*[nOrientation];
  rbM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    lbM[s] = new fullMatrix<double>((*orientationStat)[s], 9 * G);

  for(unsigned int s = 0; s < nOrientation; s++)
    rbM[s] = new fullMatrix<double>((*orientationStat)[s], 3 * G);


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
      const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
        jac->getInvertJacobian(*element[e]);

      // Reset Gauss Point Counter
      k = 0;

      // Loop on Gauss Points
      for(int g = 0; g < G; g++){
        for(unsigned int a = 0; a < 3; a++){
          for(unsigned int b = 0; b < 3; b++){
            (*lbM[s])(j, k) = 0;

            for(unsigned int i = 0; i < 3; i++)
              (*lbM[s])(j, k) +=
                (*invJac[g]->first)(i, a) *
                (*invJac[g]->first)(i, b);

            (*lbM[s])(j, k) *= fabs(invJac[g]->second);

            k++;
          }
        }
      }

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

        (*rbM[s])(j, g * 3)     = 0;
        (*rbM[s])(j, g * 3 + 1) = 0;
        (*rbM[s])(j, g * 3 + 2) = 0;

        for(unsigned int i = 0; i < 3; i++){
          (*rbM[s])(j, g * 3)     += (*invJac[g]->first)(i, 0) * fxyz(i);
          (*rbM[s])(j, g * 3 + 1) += (*invJac[g]->first)(i, 1) * fxyz(i);
          (*rbM[s])(j, g * 3 + 2) += (*invJac[g]->first)(i, 2) * fxyz(i);
        }

        (*rbM[s])(j, g * 3)     *= fabs(invJac[g]->second);
        (*rbM[s])(j, g * 3 + 1) *= fabs(invJac[g]->second);
        (*rbM[s])(j, g * 3 + 2) *= fabs(invJac[g]->second);
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}

void FormulationProjectionVector::computeA(void){
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

void FormulationProjectionVector::deleteCB(void){
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
double FormulationProjectionVector::weak(int dofI, int dofJ,
					 const GroupOfDof& god) const{
  // Init //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  // fullMatrix<double> invJac(3, 3);

  double integral = 0;
  double det;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  // MElement&      celement = const_cast<MElement&>(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
    jac->getInvertJacobian(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    // det = celement.getJacobian((*gC)(g, 0),
    //                            (*gC)(g, 1),
    //                            (*gC)(g, 2),
    //                            invJac);
    // invJac.invertInPlace();

    det = invJac[g]->second;

    phiI = Mapper::grad(eFun(dofI, g * 3),
			eFun(dofI, g * 3 + 1),
			eFun(dofI, g * 3 + 2),
			*invJac[g]->first);

    phiJ = Mapper::grad(eFun(dofJ, g * 3),
			eFun(dofJ, g * 3 + 1),
			eFun(dofJ, g * 3 + 2),
			*invJac[g]->first);

    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return integral;
}
*/
/*
double FormulationProjectionVector::rhs(int equationI,
					const GroupOfDof& god) const{
  // Init //
  fullVector<double> phi(3);
  double det;

  fullVector<double> xyz(3);
  SPoint3            pxyz;
  fullVector<double> fxyz;

  double integral = 0;
  // fullMatrix<double> invJac(3, 3);

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
    jac->getInvertJacobian(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    // Compute phi
    det = invJac[g]->second;

    // det = celement.getJacobian((*gC)(g, 0),
    //                            (*gC)(g, 1),
    //                            (*gC)(g, 2),
    //                            invJac);
    // invJac.invertInPlace();

    phi = Mapper::grad(eFun(equationI, g * 3),
		       eFun(equationI, g * 3 + 1),
		       eFun(equationI, g * 3 + 2),
		       *invJac[g]->first);

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
