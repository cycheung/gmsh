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
  WriterMsh writer;  
  
  // Get Mesh //
  Mesh msh(argv[1]);
  GroupOfElement goe = msh.getFromPhysical(7);
  //cout << msh.toString()                    << endl;
  //cout << msh.getFromPhysical(7).toString() << endl;

  /*
  // Laplace //
  FormulationLaplace laplace(goe);
  System sysLaplace(laplace);
  
  sysLaplace.fixBC(msh.getFromPhysical(6), -1);
  sysLaplace.fixBC(msh.getFromPhysical(5),  2);

  sysLaplace.assemble();
  //sysLaplace.getMatrix().print();
  //sysLaplace.getRHS().print();
  sysLaplace.solve();
  //sysLaplace.getSol().print();

  Solution solLaplace(sysLaplace, msh);
  solLaplace.write("laplace", writer);
  */
    
  // Projection //
  fullVector<double> f(2); 
  f(0) = 1; f(1) = -1; // Vector to project
  
  FormulationProjection projection(goe, f);
  System sysProj(projection);

  sysProj.assemble();
  //sysProj.getMatrix().print();
  //sysProj.getRHS().print();
  sysProj.solve();
  //sysProj.getSol().print();
  
  Solution solProj(sysProj, msh);
  solProj.write("projection", writer);
  

  // Stop Gmsh //
  GmshFinalize();
  
  return 0;
}

