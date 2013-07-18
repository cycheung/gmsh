#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "Interpolator.h"
#include "WriterMsh.h"

#include "FormulationSteadyWaveVector.h"
#include "FormulationSteadyWaveVectorSlow.h"

#include "Timer.h"
#include "SmallFem.h"

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

void compute(const Options& option){
  // Start Timer //
  Timer timer, assemble, solve;
  timer.start();

  // Writer //
  WriterMsh writer;

  // Get Domains //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement source = msh.getFromPhysical(5);
  GroupOfElement wall   = msh.getFromPhysical(6);

  // Get Parameters //
  const double puls  = atof(option.getValue("-k")[0].c_str());
  const size_t order = atoi(option.getValue("-o")[0].c_str());

  // SteadyWaveVector //
  FormulationSteadyWaveVector sWave(domain, puls * 1, order);
  System sys(sWave);

  sys.dirichlet(source, fSource);
  sys.dirichlet(wall,   fWall);

  //sys.fixCoef(source, 1);
  //sys.fixCoef(wall,   0);

  cout << "Vectorial Steady Wave (Order: " << order
       << " --- Pulsation: "               << puls
       << "): " << sys.getSize()           << endl;

  assemble.start();
  sys.assemble();
  assemble.stop();
  cout << "Assembled: " << assemble.time() << assemble.unit()
       << endl << flush;

  solve.start();
  sys.solve();
  solve.stop();
  cout << "Solved: " << solve.time() << solve.unit()
       << endl << flush;

  if(!option.getValue("-nopos").size()){
    if(option.getValue("-interp").size()){
      // Interpolated View //
      // Visu Mesh
      Mesh visuMsh(option.getValue("-interp")[0]);
      GroupOfElement visu = visuMsh.getFromPhysical(7);

      Interpolator interp(sys, visu);
      interp.write("swavev", writer);
    }

    else{
      // Adaptive View //
      writer.setValues(sys);
      writer.write("swavev");
    }
  }

  // Timer -- Finalize -- Return //
  timer.stop();

  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
