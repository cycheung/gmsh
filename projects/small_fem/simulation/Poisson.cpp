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
	      GroupOfElement& constraintDomain, 
	      Writer& writer, int order);

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer; 
  
  // Get Domains //
  Mesh msh(argv[1]);
  GroupOfElement           domain = msh.getFromPhysical(7);
  GroupOfElement constraintDomain = msh.getFromPhysical(5);

  cout << "Number of Element: " << domain.getNumber() 
       << endl << flush;

  // Compute FEM //
  unsigned int order = atoi(argv[2]);
  fPoisson(domain, 
	   constraintDomain,
	   writer, 
	   order);

  GmshFinalize();
  return 0;
}

void fPoisson(GroupOfElement& domain, 
	      GroupOfElement& constraintDomain,
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
  Solution solPoisson(sysPoisson);
  solPoisson.write("poisson", writer);

  // Adaptive View //
  writer.setValues(sysPoisson);
  writer.write("aPoisson");
}
