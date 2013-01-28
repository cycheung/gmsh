#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "Interpolator.h"
#include "WriterMsh.h"

#include "FormulationSteadyWaveVector.h"
#include "SlowFormulationSteadyWaveVector.h"

#include "Gmsh.h"

using namespace std;

double pi = 3.14159265359;

fullVector<double> fSource(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = 0;
  res(1) = 1;
  res(2) = 0;

  return res;
}

fullVector<double> fWall(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = 0;
  res(1) = 0;
  res(2) = 0;

  return res;
}

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

  // SteadyWaveVector //
  Formulation* sWave;

  if(argc == 5 && !strcmp(argv[4], "slow")){
    cout << "Slow Version" << endl;
    sWave = new SlowFormulationSteadyWaveVector(domain, puls * 1, order);
  }

  else{
    cout << "Fast Version" << endl;
    sWave = new FormulationSteadyWaveVector(domain, puls * 1, order);
  }

  System sys(*sWave);

  sys.dirichlet(source, fSource);
  sys.dirichlet(wall,   fWall);

  cout << "Vectorial Steady Wave (Order: " << order
       << " --- Pulsation: "               << puls
       << "): " << sys.getSize()           << endl;

  sys.assemble();
  sys.solve();

  if(argc == 5 && strcmp(argv[4], "slow")){
    // Interpolated View //
    // Visu Mesh
    Mesh visuMsh(argv[4]);
    GroupOfElement visu = visuMsh.getFromPhysical(7);

    Interpolator interp(sys, visu);
    interp.write("swavev", writer);
  }

  else{
    // Adaptive View //
    writer.setValues(sys);
    writer.write("swavev");
  }

  GmshFinalize();
  return 0;
}
