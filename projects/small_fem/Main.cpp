#include <iostream>
#include <string>
#include <vector>
#include "Mesh.h"
#include "fullMatrix.h"
#include "FormulationLaplace.h"
#include "FormulationProjection.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"

#include "Gmsh.h"

using namespace std;

int main(int argc, char** argv){
  // Init Gmsh //
  GmshInitialize(argc, argv);
  
  // Get Mesh //
  Mesh msh(argv[1]);

  FormulationLaplace laplace(*msh.getFromPhysical(7).at(0));

  System sysLaplace(laplace);
  sysLaplace.fixBC(*msh.getFromPhysical(6).at(0), -1);
  sysLaplace.fixBC(*msh.getFromPhysical(6).at(1), -1);
  sysLaplace.fixBC(*msh.getFromPhysical(6).at(2), -1);
  sysLaplace.fixBC(*msh.getFromPhysical(6).at(3), -1);

  sysLaplace.fixBC(*msh.getFromPhysical(5).at(0),  2);
  sysLaplace.fixBC(*msh.getFromPhysical(5).at(1),  2);
  sysLaplace.fixBC(*msh.getFromPhysical(5).at(2),  2);
  sysLaplace.fixBC(*msh.getFromPhysical(5).at(3),  2);
  
  sysLaplace.assemble();

  //sysLaplace.getMatrix().print();
  //sysLaplace.getRHS().print();

  sysLaplace.solve();

  //sysLaplace.getSol().print();

  Solution solLaplace(sysLaplace);

  WriterMsh writer;  
  solLaplace.write("laplace", writer);

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

