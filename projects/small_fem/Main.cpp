#include <iostream>

#include "Mesh.h"
#include "fullMatrix.h"
#include "FormulationLaplace.h"
#include "FormulationProjection.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"
#include "WriterVector.h"

#include "BasisTest.h"
#include "Gmsh.h"

using namespace std;

int run(int argc, char** argv);

int main(int argc, char** argv){
  int ret;

  //ret = basisTest(argc, argv);
  ret = run(argc, argv);

  return ret;
}

int run(int argc, char** argv){
  // Init Gmsh //
  GmshInitialize(argc, argv);
  WriterMsh    mWriter;  
  WriterVector vWriter;
  

  // Get Mesh //
  Mesh msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);
  //cout << msh.toString() << endl;

  // Laplace //
  FormulationLaplace laplace(domain, 1);
  System sysLaplace(laplace);

  sysLaplace.fixBC(msh.getFromPhysical(6), -1);
  sysLaplace.fixBC(msh.getFromPhysical(5),  2);

  sysLaplace.assemble();
  sysLaplace.solve();

  Solution solLaplace(sysLaplace);
  solLaplace.write("laplace", mWriter);
    
  
  // Projection //
  fullVector<double> f(3); 
  f(0) =  0.5; 
  f(1) = -1; 
  f(2) =  0; // Vector to project

  FormulationProjection projection(domain, f);
  System sysProj(projection);

  sysProj.assemble();
  sysProj.solve();

  Solution solProj(sysProj);
  solProj.write("projection", mWriter);  
  
  // Stop Gmsh //
  GmshFinalize();
      
  return 0;
}
