#include <cmath>

#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Mapper.h"

#include "Exception.h"
#include "FormulationLaplace.h"

using namespace std;

FormulationLaplace::FormulationLaplace(GroupOfElement& goe,
				       unsigned int order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Laplace formulation of order 0");

  // Function Space & Basis //
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Gaussian Quadrature Data //
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), 2 * (order - 1), *gC, *gW);

  G = gW->size(); // Nbr of Gauss points

  // PreEvaluate
  basis->preEvaluateDerivatives(*gC);

  // Fast Assembly //
  goe.orientAllElements(*basis);
  jac = new Jacobian(goe, *gC);
  jac->computeInvertJacobians();
  this->goe = &goe;

  orientationStat = &goe.getOrientationStats();
  nOrientation    = basis->getNOrientation();
  nFunction       = basis->getNFunction();
  nElement        = goe.getNumber();

  eMap = new map<const MElement*, pair<unsigned int, unsigned int> >;

  computeC();
  computeB();
  computeA();

  deleteCB();
}

FormulationLaplace::~FormulationLaplace(void){
  delete gC;
  delete gW;
  delete basis;
  delete fspace;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete aM[s];

  delete[] aM;

  delete eMap;
}
/*
double FormulationLaplace::weak(int dofI, int dofJ,
				const GroupOfDof& god) const{
  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double> invJac(3, 3);
  double integral = 0;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedDerivatives(element);

  // Loop over Integration Point //
  for(int g = 0; g < G; g++){
    double det = celement.getJacobian((*gC)(g, 0),
				      (*gC)(g, 1),
				      (*gC)(g, 2),
				      invJac);
    invJac.invertInPlace();

    phiI = Mapper::grad(eFun(dofI, g * 3),
			eFun(dofI, g * 3 + 1),
			eFun(dofI, g * 3 + 2),
			invJac);

    phiJ = Mapper::grad(eFun(dofJ, g * 3),
			eFun(dofJ, g * 3 + 1),
			eFun(dofJ, g * 3 + 2),
			invJac);

    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return integral;
}
*/

double FormulationLaplace::weak(int dofI, int dofJ,
				const GroupOfDof& god) const{

  map<const MElement*, pair<unsigned int, unsigned int> >::iterator
    index = eMap->find(&god.getGeoElement());

  return (*aM[index->second.first])
    (index->second.second, dofI * nFunction + dofJ);
}

void FormulationLaplace::computeC(void){
  unsigned int k;
  unsigned int l;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(9 * G, nFunction * nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis->getPreEvaluatedDerivatives(s);

    // Reset Gauss Point Counter
    k = 0;

    // Loop on Gauss Points
    for(int g = 0; g < G; g++){
      for(unsigned int a = 0; a < 3; a++){
        for(unsigned int b = 0; b < 3; b++){
          // Reset Function Counter
          l = 0;

          // Loop on Functions
          for(unsigned int i = 0; i < nFunction; i++){
            for(unsigned int j = 0; j < nFunction; j++){
              (*cM[s])(k, l) = (*gW)(g) * phi(i, g * 3 + a) * phi(j, g * 3 + b);
              l++;
            }
          }

          k++;
        }
      }
    }
  }
}

void FormulationLaplace::computeB(void){
  unsigned int offset = 0;
  unsigned int j;
  unsigned int k;

  // Alloc //
  bM = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    bM[s] = new fullMatrix<double>((*orientationStat)[s], 9 * G);


  // Fill //
  const vector<const MElement*>& element = goe->getAll();

  for(unsigned int s = 0; s < nOrientation; s++){
    // Reset Element Counter
    j = 0;

    // Loop on Elements
    for(unsigned int e = offset; e < offset + (*orientationStat)[s]; e++){
      // Add to eMap
      eMap->insert(pair<const MElement*, pair<unsigned int, unsigned int> >
                     (element[e], pair<unsigned int, unsigned int>(s, j)));

      // Get Jacobians
      const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
        jac->getInvertJacobian(*element[e]);

      // Reset Gauss Point Counter
      k = 0;

      // Loop on Gauss Points
      for(int g = 0; g < G; g++){
        for(unsigned int a = 0; a < 3; a++){
          for(unsigned int b = 0; b < 3; b++){
            (*bM[s])(j, k) = 0;

            for(unsigned int i = 0; i < 3; i++)
              (*bM[s])(j, k) +=
                (*invJac[g]->first)(i, a) *
                (*invJac[g]->first)(i, b);

            (*bM[s])(j, k) *= fabs(invJac[g]->second);
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

void FormulationLaplace::computeA(void){
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

void FormulationLaplace::deleteCB(void){
  for(unsigned int s = 0; s < nOrientation; s++)
    delete cM[s];

  delete[] cM;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete bM[s];

  delete[] bM;
}
