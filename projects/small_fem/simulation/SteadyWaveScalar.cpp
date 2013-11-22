#include <iostream>

#include "Mesh.h"
#include "System.h"
#include "Interpolator.h"
#include "WriterMsh.h"

#include "FormulationSteadyWaveScalar.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

double fSource(fullVector<double>& xyz){
  return 1;
}

double fWall(fullVector<double>& xyz){
  return 0;
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
  assemble.start();
  FormulationSteadyWaveScalar sWave(domain, puls * 1, order);
  System sys(sWave);

  sys.dirichlet(source, fSource);
  sys.dirichlet(wall,   fWall);

  //sys.fixCoef(source, 1);
  //sys.fixCoef(wall,   0);

  cout << "Scalar Steady Wave (Order: " << order
       << " --- Pulsation: "               << puls
       << "): " << sys.getSize()           << endl;

  sys.assemble();
  assemble.stop();
  cout << "Assembled: " << assemble.time() << assemble.unit()
       << endl << flush;

  solve.start();
  sys.solve();
  solve.stop();
  cout << "Solved: " << solve.time() << solve.unit()
       << endl << flush;

  // Write Sol //
  if(!option.getValue("-nopos").size()){
    if(option.getValue("-interp").size()){
      // Interpolated View //
      // Visu Mesh
      Mesh visuMsh(option.getValue("-interp")[0]);
      GroupOfElement visu = visuMsh.getFromPhysical(7);

      Interpolator interp(sys, visu);
      interp.write("swaves", writer);
    }

    else{
      // Adaptive View //
      writer.setValues(sys);
      writer.write("swaves");
    }
  }

  // Timer -- Finalize -- Return //
  timer.stop();

  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-nopos,-interp");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
