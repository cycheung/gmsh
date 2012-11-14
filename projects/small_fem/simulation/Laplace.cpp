#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"

#include "FormulationLaplace.h"

#include "Gmsh.h"

using namespace std;

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;
  
  // Get Mesh //
  Mesh msh(argv[1]);
 
  // Get Domain //
  GroupOfElement domain = msh.getFromPhysical(7);

  // Laplace //  
  FormulationLaplace laplace(domain, 1);
  System sysLaplace(laplace);

  sysLaplace.fixCoef(msh.getFromPhysical(6), -1);
  sysLaplace.fixCoef(msh.getFromPhysical(5),  2);

  sysLaplace.assemble();
  cout << "Laplace: " << sysLaplace.getSize() << endl;
  sysLaplace.solve();

  Solution solLaplace(sysLaplace);
  solLaplace.write("laplace", writer);

  GmshFinalize();
  return 0;
}
