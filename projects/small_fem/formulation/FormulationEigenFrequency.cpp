#include "Exception.h"
#include "fullMatrix.h"
#include "GaussIntegration.h"
#include "Polynomial.h"
#include "Mapper.h"

#include "FormulationEigenFrequency.h"

using namespace std;

// Pi  = atan(1) * 4
// Mu  = 4 * Pi * 10^-7
// Eps = 8.85418781762 * 10^−12
//const double FormulationEigenFrequency::mu  = 4 * atan(1) * 4 * 1E-7;
//const double FormulationEigenFrequency::eps = 8.85418781762E-12;

const double FormulationEigenFrequency::mu  = 1;
const double FormulationEigenFrequency::eps = 1;

FormulationEigenFrequency::FormulationEigenFrequency(const GroupOfElement& goe,
						     unsigned int order){
  // Function Space //
  fspace = new FunctionSpaceEdge(goe, order);

  // Gaussian Quadrature Data (Term One) // 
  // NB: We need to integrad a rot * rot !
  //     and order(rot f) = order(f) - 1
  gC1 = new fullMatrix<double>();
  gW1 = new fullVector<double>();

  // Gaussian Quadrature Data (Term Two) //
  // NB: We need to integrad a f * f !
  gC2 = new fullMatrix<double>();
  gW2 = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  
  // if Order == 0 --> we want Nedelec Basis of ordre *almost* one //
  if(order == 0){
    gaussIntegration::get(goe.get(0).getType(), 0, *gC1, *gW1);
    gaussIntegration::get(goe.get(0).getType(), 2, *gC2, *gW2);
  }

  else{
    gaussIntegration::get(goe.get(0).getType(), (order - 1) + (order - 1), *gC1, *gW1);
    gaussIntegration::get(goe.get(0).getType(), order + order, *gC2, *gW2);
  }

  // Nbr of Gauss points
  G1 = gW1->size(); 
  G2 = gW2->size(); 

  // PreEvaluate
  fspace->preEvaluateCurlLocalFunctions(*gC1);
  fspace->preEvaluateLocalFunctions(*gC2);
}

FormulationEigenFrequency::~FormulationEigenFrequency(void){
  delete gC1;
  delete gW1;
  delete gC2;
  delete gW2;
  delete fspace;
}

double FormulationEigenFrequency::weakA(int dofI, int dofJ,
					const GroupOfDof& god) const{
  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double> jac(3, 3);        

  double integral = 0;
  double det;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const fullMatrix<double>& eFun = 
    fspace->getEvaluatedCurlLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G1; g++){
    det = celement.getJacobian((*gC1)(g, 0), 
			       (*gC1)(g, 1), 
			       (*gC1)(g, 2), 
			       jac);

    phiI = Mapper::curl(eFun(dofI, g * 3),
			eFun(dofI, g * 3 + 1),
			eFun(dofI, g * 3 + 2), 
			jac, 1 / det);
    
    phiJ = Mapper::curl(eFun(dofJ, g * 3),
			eFun(dofJ, g * 3 + 1),
			eFun(dofJ, g * 3 + 2), 
			jac, 1 / det);

    integral += ((phiI * phiJ) / mu) * fabs(det) * (*gW1)(g);
  }

  return integral;
}

double FormulationEigenFrequency::weakB(int dofI, int dofJ,
					const GroupOfDof& god) const{
  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiJ(3);
  fullMatrix<double> invJac(3, 3);       

  double integral = 0;
  double det;

  // Get Element and Basis Functions //
  const MElement& element = god.getGeoElement();
  MElement&      celement = const_cast<MElement&>(element);
  
  const fullMatrix<double>& eFun = 
    fspace->getEvaluatedLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G2; g++){
    det = celement.getJacobian((*gC2)(g, 0), 
			       (*gC2)(g, 1), 
			       (*gC2)(g, 2), 
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

    integral += ((phiI * phiJ) * eps) * fabs(det) * (*gW2)(g);
  }

  return integral;
}
