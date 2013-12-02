#include <iostream>

#include "Mesh.h"
#include "System.h"

#include "FormulationSteadyWaveScalar.h"
#include "FormulationSteadyWaveVector.h"
#include "FormulationSteadyWaveVectorSlow.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

fullVector<double> fSourceVec(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = 0;
  res(1) = 1;
  res(2) = 0;

  return res;
}

fullVector<double> fWallVec(fullVector<double>& xyz){
  fullVector<double> res(3);

  res(0) = 0;
  res(1) = 0;
  res(2) = 0;

  return res;
}

double fSourceScal(fullVector<double>& xyz){
  return 1;
}

double fWallScal(fullVector<double>& xyz){
  return 0;
}

void compute(const Options& option){
  // Start Timer //
  Timer timer, assemble, solve;
  timer.start();

  // Get Domains //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement source = msh.getFromPhysical(5);
  GroupOfElement wall   = msh.getFromPhysical(6);

  // Get Parameters //
  const double puls  = atof(option.getValue("-k")[0].c_str());
  const size_t order = atoi(option.getValue("-o")[0].c_str());


  // Chose write formulation for Steady Wave and boundary condition //
  FormulationTyped<double>* wave = NULL;
  System*                   sys  = NULL;

  if(option.getValue("-type")[0].compare("vector") == 0){
    assemble.start();
    wave = new FormulationSteadyWaveVector(domain, puls * 1, order);
    sys  = new System(*wave);

    sys->dirichlet(source, fSourceVec);
    sys->dirichlet(wall,   fWallVec);
    cout << "Vectorial ";
  }

  else if(option.getValue("-type")[0].compare("slow") == 0){
    assemble.start();
    wave = new FormulationSteadyWaveVectorSlow(domain, puls * 1, order);
    sys  = new System(*wave);

    sys->dirichlet(source, fSourceVec);
    sys->dirichlet(wall,   fWallVec);
    cout << "Slow Vectorial ";
  }

  else if(option.getValue("-type")[0].compare("scalar") == 0){
    assemble.start();
    wave = new FormulationSteadyWaveScalar(domain, puls * 1, order);
    sys  = new System(*wave);

    sys->dirichlet(source, fSourceScal);
    sys->dirichlet(wall,   fWallScal);
    cout << "Scalar ";
  }

  else
    throw Exception("No -type given");

  cout << "Steady Wave (Order: "  << order
       << " --- Pulsation: "      << puls
       << "): " << sys->getSize() << endl;

  // Assemble and solve //
  sys->assemble();
  assemble.stop();
  cout << "Assembled: " << assemble.time() << assemble.unit()
       << endl << flush;

  solve.start();
  sys->solve();
  solve.stop();
  cout << "Solved: " << solve.time() << solve.unit()
       << endl << flush;

  // Write Sol //
  if(!option.getValue("-nopos").size()){
    FEMSolution feSol;
    sys->addSolution(feSol);
    feSol.write("swave");
  }

  // Clean //
  delete sys;
  delete wave;

  // Timer -- Finalize -- Return //
  timer.stop();

  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-nopos,-type");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
