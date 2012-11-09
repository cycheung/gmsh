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

  // Get Order //
  unsigned int order = atoi(argv[2]);
 
  // Compute //
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

  // Interpolated View //
  if(argc == 4){
    Mesh visuMesh(argv[3]);
    GroupOfElement visu = visuMesh.getFromPhysical(7);
    
    Solution sol(sysPoisson, visu);
    sol.write("iPoisson", writer);
  }

  // Adaptive View //
  writer.setValues(sysPoisson);
  writer.write("poisson");

  GmshFinalize();
  return 0;
}
