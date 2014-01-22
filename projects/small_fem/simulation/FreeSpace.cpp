#include <iostream>
#include <complex>

#include "Mesh.h"
#include "System.h"
#include "SystemHelper.h"

#include "FormulationSteadyWaveScalar.h"
#include "FormulationNeumann.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

complex<double> fSourceScal(fullVector<double>& xyz){
  return complex<double>(1, 0);
}

void compute(const Options& option){
  // Start Timer //
  Timer timer, assemble, solve;
  timer.start();

  // Get Domains //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain     = msh.getFromPhysical(7);
  GroupOfElement source     = msh.getFromPhysical(5);
  GroupOfElement freeSpace  = msh.getFromPhysical(6);

  // Get Parameters //
  const double k     = atof(option.getValue("-k")[0].c_str());
  const size_t order = atoi(option.getValue("-o")[0].c_str());

  // Formulation //
  assemble.start();
  FormulationSteadyWaveScalar<complex<double> > wave(domain, k, order);
  FormulationNeumann neumann(freeSpace, k, order);

  // System //
  System<complex<double> > sys(wave);
  SystemHelper<complex<double> >::dirichlet(sys, source, fSourceScal);

  cout << "Free Space (Order: "  << order
       << " --- Wavenumber: "    << k
       << "): " << sys.getSize() << endl;

  // Assemble and Neumann//
  sys.assemble();
  sys.addBorderTerm(neumann);
  assemble.stop();

  cout << "Assembled: " << assemble.time() << assemble.unit()
       << endl << flush;

  // Solve //
  solve.start();
  sys.solve();
  solve.stop();

  cout << "Solved: " << solve.time() << solve.unit()
       << endl << flush;

  // Write Sol //
  if(!option.getValue("-nopos").size()){
    FEMSolution<complex<double> > feSol;
    sys.getSolution(feSol);
    feSol.write("free");
  }

  // Timer -- Finalize -- Return //
  timer.stop();

  cout << "Elapsed Time: " << timer.time()
       << " s"             << endl;
}

int main(int argc, char** argv){
  // Init SmallFem //
  SmallFem::Keywords("-msh,-o,-k,-nopos");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
