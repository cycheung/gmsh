#include "Exception.h"
#include "TermProjectionGrad.h"

using namespace std;

TermProjectionGrad::
TermProjectionGrad(const GroupOfJacobian& goj,
                   const Basis& basis,
                   const fullVector<double>& integrationWeights,
                   const fullMatrix<double>& integrationPoints,
                   fullVector<double> (*f)(fullVector<double>& xyz)){
  // Basis Check //
  bFunction getFunction;

  switch(basis.getType()){
  case 0:
    getFunction = &Basis::getPreEvaluatedDerivatives;
    break;

  case 1:
    getFunction = &Basis::getPreEvaluatedFunctions;
    break;

  default:
    throw
      Exception
      ("A Grad Term must use a 1form basis, or a (gradient of) 0form basis");
  }

  // Orientations & Functions //
  orientationStat = &goj.getAllElements().getOrientationStats();
  nOrientation    = basis.getReferenceSpace().getNReferenceSpace();
  nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

  computeC(basis, getFunction, integrationWeights, cM);
  computeB(goj, basis, integrationPoints, f, bM);

  allocA(nFunction);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

TermProjectionGrad::~TermProjectionGrad(void){
}

void TermProjectionGrad::computeC(const Basis& basis,
                                  const bFunction& getFunction,
                                  const fullVector<double>& gW,
                                  fullMatrix<double>**& cM){

  const size_t nG = gW.size();
  size_t k;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(3 * nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      (basis.*getFunction)(s);

    // Loop on Gauss Points
    k = 0;

    for(size_t g = 0; g < nG; g++){
      for(size_t a = 0; a < 3; a++){
        for(size_t i = 0; i < nFunction; i++)
          (*cM[s])(k, i) = gW(g) * phi(i, k);

        k++;
      }
    }
  }
}

void TermProjectionGrad::
computeB(const GroupOfJacobian& goj,
         const Basis& basis,
         const fullMatrix<double>& gC,
         fullVector<double> (*f)(fullVector<double>& xyz),
         fullMatrix<double>**& bM){

  const size_t nG = gC.size1();
  size_t offset = 0;
  size_t j;

  fullVector<double> xyz(3);
  double             pxyz[3];
  fullVector<double> fxyz;

  // Alloc //
  bM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    bM[s] = new fullMatrix<double>((*orientationStat)[s], 3 * nG);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Loop On Element
    j = 0;

    for(size_t e = offset; e < offset + (*orientationStat)[s]; e++){
      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
        goj.getJacobian(e).getInvertJacobianMatrix();

      // Loop on Gauss Points
      for(size_t g = 0; g < nG; g++){
        // Compute f in the *physical* coordinate
        basis.getReferenceSpace().mapFromABCtoXYZ(goj.getAllElements().get(e),
                                                  gC(g, 0),
                                                  gC(g, 1),
                                                  gC(g, 2),
                                                  pxyz);
        xyz(0) = pxyz[0];
        xyz(1) = pxyz[1];
        xyz(2) = pxyz[2];

        fxyz = f(xyz);

        // Compute B
        (*bM[s])(j, g * 3)     = 0;
        (*bM[s])(j, g * 3 + 1) = 0;
        (*bM[s])(j, g * 3 + 2) = 0;

        for(size_t i = 0; i < 3; i++){
          (*bM[s])(j, g * 3)     += (*invJac[g]->first)(i, 0) * fxyz(i);
          (*bM[s])(j, g * 3 + 1) += (*invJac[g]->first)(i, 1) * fxyz(i);
          (*bM[s])(j, g * 3 + 2) += (*invJac[g]->first)(i, 2) * fxyz(i);
        }

        (*bM[s])(j, g * 3)     *= fabs(invJac[g]->second);
        (*bM[s])(j, g * 3 + 1) *= fabs(invJac[g]->second);
        (*bM[s])(j, g * 3 + 2) *= fabs(invJac[g]->second);
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}
