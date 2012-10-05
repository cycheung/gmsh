#include <iostream>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"

#include "FormulationPoisson.h"

#include "Gmsh.h"

using namespace std;

void fPoisson(Mesh& msh, Mesh& visu, Writer& writer, int order);

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer; 
  
  // Get Mesh //
  Mesh msh(argv[1]);
  Mesh visu(argv[2]);

  GroupOfElement domain = msh.getFromPhysical(9);
  cout << "Number of Element: " << domain.getNumber() << endl;

  // Compute FEM //
  unsigned int order = atoi(argv[3]);
  fPoisson(msh, visu,  writer, order);

  GmshFinalize();
  return 0;
}

void fPoisson(Mesh& msh, Mesh& visu, Writer& writer, int order){
  // FEM Solution
  GroupOfElement domain = msh.getFromPhysical(9);

  FormulationPoisson poisson(domain, order);
  System sysPoisson(poisson);
  
  sysPoisson.fixDof(msh.getFromPhysical(5), 0);
  sysPoisson.fixDof(msh.getFromPhysical(6), 0);
  sysPoisson.fixDof(msh.getFromPhysical(7), 0);
  sysPoisson.fixDof(msh.getFromPhysical(8), 0);
  
  sysPoisson.assemble();
  cout << "Poisson (" << order << "): " << sysPoisson.getSize() << endl;
  //cout << "Function Space:" << endl << poisson.fs().toString() << endl;
  sysPoisson.solve();

  GroupOfElement visuDomain = visu.getFromPhysical(9);
  Solution solPoisson(sysPoisson, visuDomain);

  solPoisson.write("poisson", writer);
}
