#include "Exception.h"
#include "TermProjectionField.h"

using namespace std;

TermProjectionField::TermProjectionField(const GroupOfJacobian& goj,
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
  nOrientation    = basis.getNOrientation();
  nFunction       = basis.getNFunction();

  // Compute //
  fullMatrix<double>** cM;
  fullMatrix<double>** bM;

  computeC(basis, integrationWeights, cM);
  computeB(goj, integrationPoints, f, bM);
  computeA(bM, cM);

  // Clean up //
  clean(bM, cM);
}

TermProjectionField::~TermProjectionField(void){
}

void TermProjectionField::computeC(const Basis& basis,
                                   const fullVector<double>& gW,
                                   fullMatrix<double>**& cM){

  const unsigned int nG = gW.size();

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(nG, nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(unsigned int g = 0; g < nG; g++)
      for(unsigned int i = 0; i < nFunction; i++)
        (*cM[s])(g, i) = gW(g) * phi(i, g);
  }
}

void TermProjectionField::computeB(const GroupOfJacobian& goj,
                                   const fullMatrix<double>& gC,
                                   double (*f)(fullVector<double>& xyz),
                                   fullMatrix<double>**& bM){

  const unsigned int nG = gC.size1();
  unsigned int offset = 0;
  unsigned int j;

  fullVector<double> xyz(3);
  SPoint3            pxyz;
  double             fxyz;

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

      for(unsigned int g = 0; g < nG; g++){
        // Compute f in the *physical* coordinate
        const_cast<MElement&>(goj.getAllElements().get(e))
          .pnt(gC(g, 0),
               gC(g, 1),
               gC(g, 2),
               pxyz);

        xyz(0) = pxyz.x();
        xyz(1) = pxyz.y();
        xyz(2) = pxyz.z();

        fxyz = f(xyz);

        (*bM[s])(j, g) = fabs(jacM[g]->second) * fxyz;
      }

      // Next Element in Orientation[s]
      j++;
    }

    // New Offset
    offset += (*orientationStat)[s];
  }
}
