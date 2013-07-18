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
  nOrientation    = basis.getReferenceSpace().getNReferenceSpace();
  nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

  computeC(basis, integrationWeights, cM);
  computeB(goj, integrationWeights.size(), bM);

  allocA(nFunction * nFunction);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

TermFieldField::~TermFieldField(void){
}

void TermFieldField::computeC(const Basis& basis,
                              const fullVector<double>& gW,
                              fullMatrix<double>**& cM){

  const size_t nG = gW.size();
  size_t l;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(nG, nFunction * nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(size_t g = 0; g < nG; g++){

      // Loop on Functions
      l = 0;

      for(size_t i = 0; i < nFunction; i++){
        for(size_t j = 0; j < nFunction; j++){
          (*cM[s])(g, l) = gW(g) * phi(i, g) * phi(j, g);
          l++;
        }
      }
    }
  }
}

void TermFieldField::computeB(const GroupOfJacobian& goj,
                              size_t nG,
                              fullMatrix<double>**& bM){
  size_t offset = 0;
  size_t j;

  // Alloc //
  bM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    bM[s] = new fullMatrix<double>((*orientationStat)[s], nG);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Loop On Element
    j = 0;

    for(size_t e = offset; e < offset + (*orientationStat)[s]; e++){
      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& jacM =
        goj.getJacobian(e).getJacobianMatrix();

      // Loop on Gauss Points
      for(size_t g = 0; g < nG; g++)
        (*bM[s])(j, g) = fabs(jacM[g]->second);;

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}
