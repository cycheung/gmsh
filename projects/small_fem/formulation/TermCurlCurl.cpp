#include "Exception.h"
#include "TermCurlCurl.h"

using namespace std;

TermCurlCurl::TermCurlCurl(const GroupOfJacobian& goj,
                           const Basis& basis,
                           const fullVector<double>& integrationWeights){
  // Basis Check //
  bFunction getFunction;

  switch(basis.getType()){
  case 1:
    getFunction = &Basis::getPreEvaluatedDerivatives;
    break;

  case 2:
    getFunction = &Basis::getPreEvaluatedFunctions;
    break;

  default:
    throw
      Exception
      ("A Curl Curl Term must use a 2form basis, or a (curl of) 1form basis");
  }

  // Orientations & Functions //
  orientationStat = &goj.getAllElements().getOrientationStats();
  nOrientation    = basis.getNOrientation();
  nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

  computeC(basis, getFunction, integrationWeights, cM);
  computeB(goj, integrationWeights.size(), bM);

  allocA(nFunction * nFunction);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

TermCurlCurl::~TermCurlCurl(void){
}

void TermCurlCurl::computeC(const Basis& basis,
                            const bFunction& getFunction,
                            const fullVector<double>& gW,
                            fullMatrix<double>**& cM){

  const unsigned int nG = gW.size();
  unsigned int k;
  unsigned int l;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(9 * nG, nFunction * nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      (basis.*getFunction)(s);

    // Loop on Gauss Points
    k = 0;

    for(unsigned int g = 0; g < nG; g++){
      for(unsigned int a = 0; a < 3; a++){
        for(unsigned int b = 0; b < 3; b++){
          // Loop on Functions
          l = 0;

          for(unsigned int i = 0; i < nFunction; i++){
            for(unsigned int j = 0; j < nFunction; j++){
              (*cM[s])(k, l) =
                gW(g) * phi(i, g * 3 + a) * phi(j, g * 3 + b);

              l++;
            }
          }

          k++;
        }
      }
    }
  }
}

void TermCurlCurl::computeB(const GroupOfJacobian& goj,
                            unsigned int nG,
                            fullMatrix<double>**& bM){
  unsigned int offset = 0;
  unsigned int j;
  unsigned int k;

  // Alloc //
  bM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    bM[s] = new fullMatrix<double>((*orientationStat)[s], 9 * nG);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Loop On Element
    j = 0;

    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& MJac =
        goj.getJacobian(e).getJacobianMatrix();

      // Loop on Gauss Points
      k = 0;

      for(unsigned int g = 0; g < nG; g++){
        for(unsigned int a = 0; a < 3; a++){
          for(unsigned int b = 0; b < 3; b++){
            (*bM[s])(j, k) = 0;

            for(unsigned int i = 0; i < 3; i++)
              (*bM[s])(j, k) +=
                (*MJac[g]->first)(i, a) * (*MJac[g]->first)(i, b);

            (*bM[s])(j, k) /= fabs(MJac[g]->second);

            k++;
          }
        }
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}
