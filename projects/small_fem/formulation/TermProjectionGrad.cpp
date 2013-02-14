#include "Exception.h"
#include "TermProjectionGrad.h"

using namespace std;

TermProjectionGrad::TermProjectionGrad(const Jacobian& jac,
                                       const Basis& basis,
                                       const fullVector<double>& integrationWeights,
                                       const fullMatrix<double>& integrationPoints,
                                       fullVector<double> (*f)(fullVector<double>& xyz)){
  // Basis Check //
  bool derivative;

  switch(basis.getType()){
  case 0:
    derivative = true;
    break;

  case 1:
    derivative = false;
    break;

  default:
    throw
      Exception
      ("A Grad Term must use a 1form basis, or a (gradient of) 0form basis");
  }

  // Function to Project //
  this->f = f;

  // Gauss Weights //
  gW = &integrationWeights;
  gC = &integrationPoints;
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

TermProjectionGrad::~TermProjectionGrad(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete aM[s];

  delete[] aM;
  delete   eMap;
}

void TermProjectionGrad::clean(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete cM[s];

  delete[] cM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete bM[s];

  delete[] bM;

  delete[] phi;
}

void TermProjectionGrad::buildEMap(void){
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

void TermProjectionGrad::computeC(void){
  unsigned int k;
  unsigned int l;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(3 * nG, nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Loop on Gauss Points
    k = 0;

    for(unsigned int g = 0; g < nG; g++){
      for(unsigned int a = 0; a < 3; a++){
        // Loop on Functions
        l = 0;

        for(unsigned int i = 0; i < nFunction; i++){
          (*cM[s])(k, l) =
            (*gW)(g) *
            (*phi[s])(i, g * 3 + a);

          l++;
        }


        k++;
      }
    }
  }
}

void TermProjectionGrad::computeB(void){
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
  const vector<pair<const MElement*, ElementData> >&
    element = jac->getAllElements().getAll();

  for(unsigned int s = 0; s < nOrientation; s++){
    // Loop On Element
    j = 0;

    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
        jac->getInvertJacobian(*element[e].first);

      // Loop on Gauss Points
      for(unsigned int g = 0; g < nG; g++){
        const_cast<MElement*>(element[e].first)
          ->pnt((*gC)(g, 0),
                (*gC)(g, 1),
                (*gC)(g, 2),
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


void TermProjectionGrad::computeA(void){
  // Alloc //
  aM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    aM[s] = new fullMatrix<double>((*orientationStat)[s], nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++)
    // GEMM doesn't like matrices with 0 Elements
    if((*orientationStat)[s])
      aM[s]->gemm(*bM[s], *cM[s]);
}
