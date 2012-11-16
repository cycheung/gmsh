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
  
  // Get Domains //
  Mesh msh(argv[1]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement source = msh.getFromPhysical(5);
  GroupOfElement wall   = msh.getFromPhysical(6);

  // Get Parameters //
  const double       puls  = atof(argv[2]);
  const unsigned int order = atoi(argv[3]);

  // SteadyWave //  
  FormulationSteadyWave sWave(domain, puls * 1, order);
  System sys(sWave);

  sys.fixCoef(source, 1);
  sys.fixCoef(wall, 0);

  cout << "Steady Wave (Order: " << order 
       << " --- Pulsation: "     << puls
       << "): " << sys.getSize() << endl;

  sys.assemble();
  sys.solve();

  if(argc == 5){
    // Interpolated View //
    // Visu Mesh
    Mesh visuMsh(argv[4]);
    GroupOfElement visu = visuMsh.getFromPhysical(7);

    Solution sol(sys, visu);
    sol.write("swave", writer);
  }

  else{
    // Adaptive View //
    writer.setValues(sys);
    writer.write("swave");
  }

  GmshFinalize();
  return 0;
}
