#include <iostream>
#include <sstream>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"

#include "Interpolator.h"
#include "WriterMsh.h"

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

  sysPoisson.fixCoef(constraintDomain, 0);

  cout << "Poisson -- Order " << order
       << ": " << sysPoisson.getSize()
       << endl << flush;

  sysPoisson.assemble();
  sysPoisson.solve();

  if(argc == 4){
  // Interpolated View //
    // Visu Mesh
    Mesh visuMesh(argv[3]);
    GroupOfElement visu = visuMesh.getFromPhysical(7);

    Interpolator interp(sysPoisson, visu);
    interp.write("poisson", writer);
  }

  else{
    // Adaptive View //
    writer.setValues(sysPoisson);
    writer.write("poisson");
  }

  GmshFinalize();
  return 0;
}
