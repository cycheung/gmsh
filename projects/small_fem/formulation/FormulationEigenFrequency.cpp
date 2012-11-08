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
  // Gaussian Quadrature Data (Term One) // 
  // NB: We need to integrad a rot * rot !
  //     and order(rot f) = order(f) - 1
  gC1 = new fullMatrix<double>();
  gW1 = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), (order - 1) + (order - 1) + 2 , *gC1, *gW1);

  G1 = gW1->size(); // Nbr of Gauss points

  // Gaussian Quadrature Data (Term Two) //
  // NB: We need to integrad a f * f !
  gC2 = new fullMatrix<double>();
  gW2 = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get(goe.get(0).getType(), order + order + 2, *gC2, *gW2);

  G2 = gW2->size(); // Nbr of Gauss points


  // Function Space //
  fspace = new FunctionSpaceEdge(goe, order);
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
  
  const vector<const vector<Polynomial>*> fun = 
    fspace->getCurlLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G1; g++){
    det = celement.getJacobian((*gC1)(g, 0), 
			       (*gC1)(g, 1), 
			       (*gC1)(g, 2), 
			       jac);

    phiI = Mapper::curl(Polynomial::at(*fun[dofI], 
				       (*gC1)(g, 0), 
				       (*gC1)(g, 1),
				       (*gC1)(g, 2)),
			jac, 1 / det);

    phiJ = Mapper::curl(Polynomial::at(*fun[dofJ], 
				       (*gC1)(g, 0), 
				       (*gC1)(g, 1), 
				       (*gC1)(g, 2)),
			jac, 1 / det);
    
    integral += phiI * phiJ * fabs(det) * (*gW1)(g) / mu;
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
  
  const vector<const vector<Polynomial>*> fun = 
    fspace->getLocalFunctions(element);

  // Loop over Integration Point //
  for(int g = 0; g < G2; g++){
    det = celement.getJacobian((*gC2)(g, 0), 
			       (*gC2)(g, 1), 
			       (*gC2)(g, 2), 
			       invJac);
    invJac.invertInPlace();

    phiI = Mapper::grad(Polynomial::at(*fun[dofI],
				       (*gC2)(g, 0), 
				       (*gC2)(g, 1),
				       (*gC2)(g, 2)),
			invJac);
    
    phiJ = Mapper::grad(Polynomial::at(*fun[dofJ],
				       (*gC2)(g, 0), 
				       (*gC2)(g, 1),
				       (*gC2)(g, 2)),
			invJac);

    integral += phiI * phiJ * fabs(det) * (*gW2)(g) * eps;
  }

  return integral;
}
