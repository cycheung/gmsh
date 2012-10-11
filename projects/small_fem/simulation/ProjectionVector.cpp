#include <iostream>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"

#include "FormulationProjectionVector.h"

#include "Gmsh.h"

using namespace std;

// Vector to Project //
fullVector<double> f(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) =  0.5;
  res(1) = -0;
  res(2) =  0;

  return res;
}

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;  
  
  // Get Mesh //
  Mesh msh(argv[1]);
  
  // Get Domain //
  GroupOfElement domain = msh.getFromPhysical(7);

  // Projection //
  FormulationProjectionVector projection(domain, f, atoi(argv[2]));
  System sysProj(projection);

  sysProj.assemble();
  cout << "Projection: " << sysProj.getSize() << endl;
  sysProj.solve();

  Solution solProj(sysProj);
  solProj.write("projection", writer);

  GmshFinalize();
  return 0;
}
