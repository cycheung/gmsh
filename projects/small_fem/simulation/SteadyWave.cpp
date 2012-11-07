#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "Solution.h"
#include "WriterMsh.h"

#include "FormulationSteadyWave.h"

#include "Gmsh.h"

using namespace std;

double pi = 3.14159265359;

int main(int argc, char** argv){
  GmshInitialize(argc, argv);

  // Writer //
  WriterMsh writer;
  
  // Get Meshes //
  Mesh msh(argv[1]);
  //Mesh visuMsh(argv[2]);

  // Get Domains //
  GroupOfElement domain = msh.getFromPhysical(7);
  //GroupOfElement visu   = visuMsh.getFromPhysical(7);

  // Get Parameters //
  const double       puls  = atof(argv[2]);
  const unsigned int order = atoi(argv[3]);

  // SteadyWave //  
  FormulationSteadyWave sWave(domain, puls * 1, order);
  System sys(sWave);

  sys.fixDof(msh.getFromPhysical(5), 1);
  sys.fixDof(msh.getFromPhysical(6), 0);

  cout << "Steady Wave (Order: " << order 
       << " --- Pulsation: "     << puls
       << "): " << sys.getSize() << endl;

  sys.assemble();
  sys.solve();

  Solution sol(sys);
  sol.write("swave", writer);

  GmshFinalize();
  return 0;
}
