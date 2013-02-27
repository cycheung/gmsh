#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "Interpolator.h"
#include "WriterMsh.h"

#include "FormulationSteadyWaveVector.h"

#include "SystemInstrumented.h"
#include "Timer.h"
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
  // Start Timer //
  Timer preCpt;
  Timer timer;

  timer.start();

  // Init Gmsh //
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
  preCpt.start();
  FormulationSteadyWaveVector sWave(domain, puls * 1, order);
  preCpt.stop();

  //SystemInstrumented sys(sWave);
  System sys(sWave);

  sys.dirichlet(source, fSource);
  sys.dirichlet(wall,   fWall);

  cout << "Vectorial Steady Wave (Order: " << order
       << " --- Pulsation: "               << puls
       << "): " << sys.getSize()           << endl;

  sys.assemble();
  sys.solve();
  /*
  cout << "Precomputing: "  << preCpt.time()   << endl;
  cout << "PreAlloc: "      << sys.preAlloc    << endl;
  cout << "Dof Lookup: "    << sys.dofLookTime << endl;
  cout << "Getting Terms: " << sys.totLHSTime    + sys.totRHSTime    << endl;
  cout << "Adding Terms: "  << sys.totAddLHSTime + sys.totAddRHSTime << endl;
  */

  if(argc == 5){
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


  // Timer -- Finalize -- Return //
  GmshFinalize();
  timer.stop();

  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;

  return 0;
}
