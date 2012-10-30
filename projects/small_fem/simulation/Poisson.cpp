#include <iostream>
#include <sstream>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"

#include "Solution.h"
#include "WriterMsh.h"
#include "Integrator.h"

#include "FormulationPoisson.h"

#include "Gmsh.h"

using namespace std;

void fPoisson(GroupOfElement& domain, 
	      GroupOfElement& visuDomain, 
	      GroupOfElement& constraintDomain, 
	      Writer& writer, int order);

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer; 
  
  // Get Mesh //
  Mesh msh(argv[1]);
  Mesh visu(argv[2]);

  GroupOfElement           domain = msh.getFromPhysical(7);
  GroupOfElement constraintDomain = msh.getFromPhysical(5);
  GroupOfElement       visuDomain = visu.getFromPhysical(7);

  cout << "Number of Element: " << domain.getNumber() 
       << endl << flush;

  // Compute FEM //
  unsigned int order = atoi(argv[3]);
  fPoisson(domain, 
	   constraintDomain,
	   visuDomain, 
	   writer, 
	   order);

  GmshFinalize();
  return 0;
}

void fPoisson(GroupOfElement& domain, 
	      GroupOfElement& constraintDomain,
	      GroupOfElement& visuDomain, 
	      Writer& writer, 
	      int order){

  // FEM Solution
  FormulationPoisson poisson(domain, order);
  System sysPoisson(poisson);

  cout << "Poisson -- Order " << order 
       << ": " << sysPoisson.getSize() 
       << endl << flush;
  
  sysPoisson.fixDof(constraintDomain, 0);
  sysPoisson.assemble();
  sysPoisson.solve();

  // Integrate Solution //
  Integrator integrator(sysPoisson);
  cout << "Integrated Solution: "
       << integrator.integrate()
       << endl;

  // Write Solution //
  Solution solPoisson(sysPoisson, visuDomain);
  solPoisson.write("poisson", writer);
}
