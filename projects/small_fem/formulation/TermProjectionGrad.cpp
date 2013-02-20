#include "Exception.h"
#include "TermProjectionGrad.h"

using namespace std;

TermProjectionGrad::TermProjectionGrad(const GroupOfJacobian& goj,
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
  nOrientation    = basis.getNOrientation();
  nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

  computeC(basis, getFunction, integrationWeights, cM);
  computeB(goj, integrationPoints, f, bM);
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

  const unsigned int nG = gW.size();
  unsigned int k;
  unsigned int l;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(3 * nG, nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      (basis.*getFunction)(s);

    // Loop on Gauss Points
    k = 0;

    for(unsigned int g = 0; g < nG; g++){
      for(unsigned int a = 0; a < 3; a++){
        // Loop on Functions
        l = 0;

        for(unsigned int i = 0; i < nFunction; i++){
          (*cM[s])(k, l) = gW(g) * phi(i, g * 3 + a);

          l++;
        }


        k++;
      }
    }
  }
}

void TermProjectionGrad::computeB(const GroupOfJacobian& goj,
                                  const fullMatrix<double>& gC,
                                  fullVector<double> (*f)(fullVector<double>& xyz),
                                  fullMatrix<double>**& bM){

  const unsigned int nG = gC.size1();
  unsigned int offset = 0;
  unsigned int j;

  fullVector<double> xyz(3);
  SPoint3            pxyz;
  fullVector<double> fxyz;

  // Alloc //
  bM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    bM[s] = new fullMatrix<double>((*orientationStat)[s], 3 * nG);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Loop On Element
    j = 0;

    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
        goj.getJacobian(e).getInvertJacobianMatrix();

      // Loop on Gauss Points
      for(unsigned int g = 0; g < nG; g++){
        const_cast<MElement&>(goj.getAllElements().get(e))
          .pnt(gC(g, 0),
               gC(g, 1),
               gC(g, 2),
               pxyz);

        xyz(0) = pxyz.x();
        xyz(1) = pxyz.y();
        xyz(2) = pxyz.z();

        fxyz = f(xyz);

        (*bM[s])(j, g * 3)     = 0;
        (*bM[s])(j, g * 3 + 1) = 0;
        (*bM[s])(j, g * 3 + 2) = 0;

        for(unsigned int i = 0; i < 3; i++){
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
