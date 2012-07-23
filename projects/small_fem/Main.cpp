#include <iostream>
#include <string>
#include <vector>
#include "Mesh.h"
#include "fullMatrix.h"
#include "FormulationLaplace.h"
#include "FormulationProjection.h"
#include "System.h"
#include "Solution.h"

#include "Gmsh.h"

using namespace std;

int main(int argc, char** argv){
  // Init Gmsh //
  GmshInitialize(argc, argv);
  
  // Get Mesh //
  Mesh msh(argv[1]);

  FormulationLaplace laplace;

  System sysLaplace(msh.getGroup(2), laplace);
  sysLaplace.assemble();

  sysLaplace.getMatrix().print();

  // Stop Gmsh //
  GmshFinalize();

  
  /*  
  // Laplace //
  FormulationLaplace laplace;
  System sysLaplace(msh.getAllNodeElements(), laplace);

  sysLaplace.assemble();

  sysLaplace.fixBC(5, -2);
  sysLaplace.fixBC(6,  1);

  sysLaplace.solve();

  Solution solLaplace(msh, laplace);
  solLaplace.write("laplace.pos", "laplace");

    
  // Projection //
  fullVector<double> f(2); 
  f(0) = -1; f(1) = 1; // Vector to project
  
  FormulationProjection projection(f);
  System sysProj(msh.getAllEdgeElements(), projection);

  sysProj.solve();
  
  Solution solProj(msh, projection);
  solProj.write("projection.pos", "projection");
  */
  
  return 0;
}

