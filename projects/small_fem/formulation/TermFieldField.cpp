#include "Exception.h"
#include "TermFieldField.h"

using namespace std;

TermFieldField::TermFieldField(const Jacobian& jac,
                               const Basis& basis,
                               const fullVector<double>& integrationWeights){

  // Basis Check //
  if(basis.getType() != 0)
    throw
      Exception
      ("A Field Field Term must use a 0form basis");

  // Gauss Points //
  gW = &integrationWeights;
  nG = gW->size();

  // Basis & Orientations //
  this->basis     = &basis;
  nOrientation    = basis.getNOrientation();
  nFunction       = basis.getNFunction();
  orientationStat = &jac.getAllElements().getOrientationStats();

  // Compute Jacobians //
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

TermFieldField::~TermFieldField(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete aM[s];

  delete[] aM;
  delete   eMap;
}

void TermFieldField::clean(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete cM[s];

  delete[] cM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete bM[s];

  delete[] bM;
}

void TermFieldField::buildEMap(void){
  const vector<const MElement*>& element = jac->getAllElements().getAll();

  eMap = new map<const MElement*, pair<unsigned int, unsigned int> >;

  unsigned int offset = 0;
  unsigned int j;

  for(unsigned int s = 0; s < nOrientation; s++){
    j = 0;

    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      eMap->insert(pair<const MElement*, pair<unsigned int, unsigned int> >
                       (element[e], pair<unsigned int, unsigned int>(s, j)));
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}

void TermFieldField::computeC(void){
  unsigned int l;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(nG, nFunction * nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis->getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(unsigned int g = 0; g < nG; g++){

      // Loop on Functions
      l = 0;

      for(unsigned int i = 0; i < nFunction; i++){
        for(unsigned int j = 0; j < nFunction; j++){
          (*cM[s])(g, l) = (*gW)(g) *phi(i, g) * phi(j, g);
          l++;
        }
      }
    }
  }
}

void TermFieldField::computeB(void){
  unsigned int offset = 0;
  unsigned int j;

  // Alloc //
  bM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    bM[s] = new fullMatrix<double>((*orientationStat)[s], nG);

  // Fill //
  const vector<const MElement*>& element = jac->getAllElements().getAll();

  for(unsigned int s = 0; s < nOrientation; s++){
    // Loop On Element
    j = 0;

    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& jacM =
        jac->getJacobian(*element[e]);

      for(unsigned int g = 0; g < nG; g++)
        (*bM[s])(j, g) = fabs(jacM[g]->second);;

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}

void TermFieldField::computeA(void){
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
