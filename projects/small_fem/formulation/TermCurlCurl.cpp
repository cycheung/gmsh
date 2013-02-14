#include "Exception.h"
#include "TermCurlCurl.h"

using namespace std;

TermCurlCurl::TermCurlCurl(const Jacobian& jac,
                           const Basis& basis,
                           const fullVector<double>& integrationWeights){

  // Basis Check //
  bool derivative;

  switch(basis.getType()){
  case 1:
    derivative = true;
    break;

  case 2:
    derivative = false;
    break;

  default:
    throw
      Exception
      ("A Curl Curl Term must use a 2form basis, or a (curl of) 1form basis");
  }

  // Gauss Weights //
  gW = &integrationWeights;
  nG = gW->size();

  // Basis & Orientations //
  nOrientation = basis.getNOrientation();
  nFunction    = basis.getNFunction();

  orientationStat = &jac.getAllElements().getOrientationStats();
  phi             = new const fullMatrix<double>*[nOrientation];

  if(derivative)
    for(unsigned int s = 0; s < nOrientation; s++)
      phi[s] = &basis.getPreEvaluatedDerivatives(s);

  else
    for(unsigned int s = 0; s < nOrientation; s++)
      phi[s] = &basis.getPreEvaluatedFunctions(s);

  // Jacobians //
  this->jac = &jac;

  // Element Map //
  buildEMap();

  // Compute //
  computeC();
  computeB();
  computeA();

  // Clean up //
  clean();
}

TermCurlCurl::~TermCurlCurl(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete aM[s];

  delete[] aM;
  delete   eMap;
}

void TermCurlCurl::clean(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete cM[s];

  delete[] cM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete bM[s];

  delete[] bM;

  delete[] phi;
}

void TermCurlCurl::buildEMap(void){
  const vector<pair<const MElement*, ElementData> >&
    element = jac->getAllElements().getAll();

  eMap = new map<const MElement*, pair<unsigned int, unsigned int> >;

  unsigned int offset = 0;
  unsigned int j;

  for(unsigned int s = 0; s < nOrientation; s++){
    j = 0;

    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      eMap->insert(pair<const MElement*, pair<unsigned int, unsigned int> >
                       (element[e].first, pair<unsigned int, unsigned int>(s, j)));
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}

void TermCurlCurl::computeC(void){
  unsigned int k;
  unsigned int l;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(9 * nG, nFunction * nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
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
                (*gW)(g) *
                (*phi[s])(i, g * 3 + a) *
                (*phi[s])(j, g * 3 + b);

              l++;
            }
          }

          k++;
        }
      }
    }
  }
}

void TermCurlCurl::computeB(void){
  unsigned int offset = 0;
  unsigned int j;
  unsigned int k;

  // Alloc //
  bM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    bM[s] = new fullMatrix<double>((*orientationStat)[s], 9 * nG);

  // Fill //
  const vector<pair<const MElement*, ElementData> >&
    element = jac->getAllElements().getAll();

  for(unsigned int s = 0; s < nOrientation; s++){
    // Loop On Element
    j = 0;

    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& MJac =
        jac->getJacobian(*element[e].first);

      // Loop on Gauss Points
      k = 0;

      for(unsigned int g = 0; g < nG; g++){
        for(unsigned int a = 0; a < 3; a++){
          for(unsigned int b = 0; b < 3; b++){
            (*bM[s])(j, k) = 0;

            for(unsigned int i = 0; i < 3; i++)
              (*bM[s])(j, k) +=
                (*MJac[g]->first)(i, a) *
                (*MJac[g]->first)(i, b);

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

void TermCurlCurl::computeA(void){
  // Alloc //
  aM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    aM[s] = new fullMatrix<double>((*orientationStat)[s], nFunction * nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++)
    // GEMM doesn't like matrices with 0 Elements
    if((*orientationStat)[s])
      aM[s]->gemm(*bM[s], *cM[s]);
}
