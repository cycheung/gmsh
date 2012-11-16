#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"

#include "FormulationLaplace.h"

#include "Gmsh.h"

using namespace std;

double f1(fullVector<double>& xyz){
  return -1;
}

double f2(fullVector<double>& xyz){
  return 2;
}

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;
  
  // Get Mesh //
  Mesh msh(argv[1]);

  // Get Order //
  unsigned int order = atoi(argv[2]);
 
  // Get Domain //
  GroupOfElement domain = msh.getFromPhysical(7);

  // Laplace //  
  FormulationLaplace laplace(domain, order);
  System sysLaplace(laplace);

  cout << "Laplace (" << order << "): " 
       << sysLaplace.getSize() << endl;

  sysLaplace.dirichlet(msh.getFromPhysical(6), f1);
  sysLaplace.dirichlet(msh.getFromPhysical(5), f2);

  sysLaplace.assemble();
  sysLaplace.solve();

  if(argc == 4){
    // Interpolated View //
    // Visu Mesh
    Mesh visuMesh(argv[3]);
    GroupOfElement visu = visuMesh.getFromPhysical(7);
    
    Solution sol(sysLaplace, visu);
    sol.write("laplace", writer);
  }

  else{
    // Adaptive View //
    writer.setValues(sysLaplace);
    writer.write("laplace");
  }

  GmshFinalize();
  return 0;
}
