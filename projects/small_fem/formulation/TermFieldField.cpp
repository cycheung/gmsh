#include "Exception.h"
#include "TermFieldField.h"

using namespace std;

TermFieldField::TermFieldField(const GroupOfJacobian& goj,
                               const Basis& basis,
                               const fullVector<double>& integrationWeights){
  // Basis Check //
  if(basis.getType() != 0)
    throw
      Exception
      ("A Field Field Term must use a 0form basis");

  // Orientations & Functions //
  orientationStat = &goj.getAllElements().getOrientationStats();
  nOrientation    = basis.getNOrientation();
  nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

  computeC(basis, integrationWeights, cM);
  computeB(goj, integrationWeights.size(), bM);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

TermFieldField::~TermFieldField(void){
}

void TermFieldField::computeC(const Basis& basis,
                              const fullVector<double>& gW,
                              fullMatrix<double>**& cM){

  const unsigned int nG = gW.size();
  unsigned int l;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(nG, nFunction * nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(unsigned int g = 0; g < nG; g++){

      // Loop on Functions
      l = 0;

      for(unsigned int i = 0; i < nFunction; i++){
        for(unsigned int j = 0; j < nFunction; j++){
          (*cM[s])(g, l) = gW(g) * phi(i, g) * phi(j, g);
          l++;
        }
      }
    }
  }
}

void TermFieldField::computeB(const GroupOfJacobian& goj,
                              unsigned int nG,
                              fullMatrix<double>**& bM){
  unsigned int offset = 0;
  unsigned int j;

  // Alloc //
  bM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    bM[s] = new fullMatrix<double>((*orientationStat)[s], nG);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Loop On Element
    j = 0;

    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& jacM =
        goj.getJacobian(e).getJacobianMatrix();

      // Loop on Gauss Points
      for(unsigned int g = 0; g < nG; g++)
        (*bM[s])(j, g) = fabs(jacM[g]->second);;

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}
