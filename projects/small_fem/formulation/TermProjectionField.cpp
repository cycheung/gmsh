#include "Exception.h"
#include "TermProjectionField.h"

using namespace std;

TermProjectionField::
TermProjectionField(const GroupOfJacobian& goj,
                    const Basis& basis,
                    const fullVector<double>& integrationWeights,
                    const fullMatrix<double>& integrationPoints,
                    double (*f)(fullVector<double>& xyz)){
  // Basis Check //
  if(basis.getType() != 0)
    throw
      Exception
      ("A Field Term must use a 0form basis");

  // Orientations & Function //
  orientationStat = &goj.getAllElements().getOrientationStats();
  nOrientation    = basis.getReferenceSpace().getNReferenceSpace();
  nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

  computeC(basis, integrationWeights, cM);
  computeB(goj, basis, integrationPoints, f, bM);

  allocA(nFunction);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

TermProjectionField::~TermProjectionField(void){
}

void TermProjectionField::computeC(const Basis& basis,
                                   const fullVector<double>& gW,
                                   fullMatrix<double>**& cM){

  const size_t nG = gW.size();

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(size_t g = 0; g < nG; g++)
      for(size_t i = 0; i < nFunction; i++)
        (*cM[s])(g, i) = gW(g) * phi(i, g);
  }
}

void TermProjectionField::computeB(const GroupOfJacobian& goj,
                                   const Basis& basis,
                                   const fullMatrix<double>& gC,
                                   double (*f)(fullVector<double>& xyz),
                                   fullMatrix<double>**& bM){

  const size_t nG = gC.size1();
  size_t offset = 0;
  size_t j;

  fullVector<double> xyz(3);
  double             pxyz[3];
  double             fxyz;

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
        (*bM[s])(j, g) = fabs(jacM[g]->second) * fxyz;
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}
