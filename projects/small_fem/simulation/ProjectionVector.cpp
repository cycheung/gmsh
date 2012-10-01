#include <iostream>

#include "Mesh.h"
#include "fullMatrix.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"

#include "FormulationProjectionVector.h"

using namespace std;

int main(int argc, char** argv){
  // Writer //
  WriterMsh writer;  
  
  // Get Mesh //
  Mesh msh(argv[1]);
  
  // Get Domain //
  GroupOfElement domain = msh.getFromPhysical(7);

  // Vector to Prject //
  fullVector<double> f(3); 
  f(0) =  0.5; 
  f(1) = -1; 
  f(2) =  0; // Vector to project

  // Projection //
  FormulationProjectionVector projection(domain, f);
  System sysProj(projection);

  sysProj.assemble();
  cout << "Projection: " << sysProj.getSize() << endl;
  sysProj.solve();

  Solution solProj(sysProj);
  solProj.write("projection", writer);
}
