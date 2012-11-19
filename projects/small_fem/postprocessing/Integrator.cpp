#include "GroupOfElement.h"
#include "GroupOfDof.h"
#include "Exception.h"

#include "FunctionSpaceNode.h"

#include "Integrator.h"

using namespace std;

Integrator::Integrator(const System& system){
  // Get FEM Solution //
  sol  = &(system.getSol());
  dofM = &(system.getDofManager());
  fs   = &(system.getFormulation().fs());
  
  // Get Mesh //
  const GroupOfElement& domain = fs->getSupport();
  
  element  = &(domain.getAll());  
  nElement = domain.getNumber();

  // Gaussian Quadrature Data // 
  gC = new fullMatrix<double>();
  gW = new fullVector<double>();

  // Look for 1st element to get element type
  // (We suppose only one type of Mesh !!)
  gaussIntegration::get((*element)[0]->getType(), fs->getOrder(), *gC, *gW);
}

Integrator::~Integrator(void){
  delete gC;
  delete gW;
}

double Integrator::integrate(void) const{
  // Select integration with respect //
  // To Function Space Type          //

  const int fsType = fs->getType();

  switch(fsType){
  case 0: return integrate(zero);
  case 1: return integrate(one);
  case 2: return integrate(two);
  case 3: return integrate(three);
    
  default: 
    throw Exception("Integrator: Unknown Function Space (type: %d)", fsType);
  }
}

double Integrator::integrate(double (*law)(MElement& element,
					   const FunctionSpace& fs, 
					   vector<double>& coef,
					   fullMatrix<double>& gC,
					   fullVector<double>& gW)) const{
  // Init //
  double integral = 0;
  
  // Interate on All Elements //
  for(unsigned int i = 0; i < nElement; i++){
    // Get GoD 
    const GroupOfDof& god = 
      fs->getGoDFromElement(*(*element)[i]);

    // Get Dof
    const vector<const Dof*>& dof  = god.getAll();
    const unsigned int        size = dof.size();
      
    // Get Coef
    vector<double> coef(size);
    for(unsigned int j = 0; j < size; j++)
      // Look in Solution
      coef[j] = 
	(*sol)(dofM->getGlobalId(*dof[j])); 

    // Const Cast For compatibility
    MElement& celement = 
      const_cast<MElement&>(*(*element)[i]);
    
    // Integrate
    integral += law(celement, *fs, coef, *gC, *gW);
  }

  // Return //
  return integral;
}

double Integrator::zero(MElement& element, 
			const FunctionSpace& fs, 
			vector<double>& coef,
			fullMatrix<double>& gC,
			fullVector<double>& gW){
  
  // If we are here fs is a *Nodal* Space //
  const FunctionSpaceNode& fNode = 
    static_cast<const FunctionSpaceNode&>(fs);

  // Init Some Stuff //
  fullMatrix<double> jac(3, 3);        
  double integral = 0;
  double det;
  double phi;
  
  // Get Local Functions //
  const vector<const Polynomial*> fun = 
    fNode.getLocalFunctions(element);

  // Get Sizes //
  const unsigned int F = fun.size();
  const unsigned int G = gW.size();
  
  // Quick Check // 
  if(F != coef.size())
    throw Exception
      ("Integrator: %d Functions for %d Coeficients !",
       F, coef.size());

  // Loop over Functions //  
  for(unsigned int f = 0; f < F; f++){
    // Loop over Integration Point
    for(unsigned int g = 0; g < G; g++){
      det = element.getJacobian(gC(g, 0), 
				gC(g, 1), 
				gC(g, 2), 
				jac);
      
      phi = fun[f]->at(gC(g, 0), 
		       gC(g, 1), 
		       gC(g, 2));
      
      integral += 
	coef[f] * phi * fabs(det) * gW(g);
    }
  }

  // Return
  return integral;
}

double Integrator::one(MElement& element, 
		       const FunctionSpace& fs, 
		       vector<double>& coef,
		       fullMatrix<double>& gC,
		       fullVector<double>& gW){

  throw Exception("Integrator: One Form not Implemented");
}

double Integrator::two(MElement& element, 
		       const FunctionSpace& fs, 
		       vector<double>& coef,
		       fullMatrix<double>& gC,
		       fullVector<double>& gW){

  throw Exception("Integrator: Two Form not Implemented");
}

double Integrator::three(MElement& element, 
			 const FunctionSpace& fs, 
			 vector<double>& coef,
			 fullMatrix<double>& gC,
			 fullVector<double>& gW){

  throw Exception("Integrator: Three Form not Implemented");
}
