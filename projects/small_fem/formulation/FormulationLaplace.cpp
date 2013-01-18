#include <cmath>

#include "BasisGenerator.h"
#include "GaussIntegration.h"
#include "Mapper.h"

#include "Exception.h"
#include "FormulationLaplace.h"

using namespace std;

FormulationLaplace::FormulationLaplace(const GroupOfElement& goe,
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
  nOrientation = basis->getNOrientation();
  nFunction    = basis->getNFunction();

  computeC();
}

FormulationLaplace::~FormulationLaplace(void){
  delete gC;
  delete gW;
  delete basis;
  delete fspace;

  for(unsigned int s = 0; s < nOrientation; s++)
    delete c[s];

  delete[] c;
}

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

void FormulationLaplace::computeC(void){
  unsigned int k;
  unsigned int l;

  // Alloc //
  c = new fullMatrix<double>*[nOrientation];

  for(unsigned int s = 0; s < nOrientation; s++)
    c[s] = new fullMatrix<double>(9 * G, nFunction * nFunction);

  // Fill //
  for(unsigned int s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis->getPreEvaluatedDerivatives(s);

    // Reset counter
    k = 0;
    l = 0;

    // Loop on Gauss Points
    for(int g = 0; g < G; g++){
      for(unsigned int a = 0; a < 3; a++){
        for(unsigned int b = 0; b < 3; b++){
          l = 0;

          // Loop on Functions
          for(unsigned int i = 0; i < nFunction; i++){
            for(unsigned int j = 0; j < nFunction; j++){
              (*c[s])(k, l) = (*gW)(g) * phi(i, g + a) * phi(j, g + b);
              l++;
            }
          }

          k++;
        }
      }
    }
  }
}
