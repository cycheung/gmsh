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
/*
fullVector<double> f(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) =  0.5;
  res(1) = -0;
  res(2) =  0;

  return res;
}
*/

fullVector<double> f(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = sin(10 * xyz(0));
  res(1) = sin(10 * xyz(1));
  res(2) = sin(10 * xyz(2));

  return res;
}


int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;  
  
  // Get Mesh //
  Mesh     msh(argv[1]);
  Mesh visuMsh(argv[2]);

  // Get Domain //
  GroupOfElement domain =     msh.getFromPhysical(7);
  GroupOfElement visu   = visuMsh.getFromPhysical(7);  

  // Get Order //
  const unsigned int order = atoi(argv[3]);

  // Projection //
  FormulationProjectionVector projection(domain, f, order);
  System sysProj(projection);

  sysProj.assemble();
  cout << "Projection: " << sysProj.getSize() << endl;
  sysProj.solve();

  Solution solProj(sysProj, visu);
  solProj.write("projection", writer);

  GmshFinalize();
  return 0;
}
