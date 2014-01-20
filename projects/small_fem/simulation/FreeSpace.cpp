#include <iostream>
#include <complex>

#include "Mesh.h"
#include "System.h"
#include "SystemHelper.h"

#include "FormulationSteadyWaveScalar.h"
#include "FormulationNeumann.h"

#include "Timer.h"
#include "SmallFem.h"

#include "NodeSolution.h"
#include "Interpolator.h"

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
  //GroupOfElement outerSpace = msh.getFromPhysical(4);

  // Get Parameters //
  const double k     = atof(option.getValue("-k")[0].c_str());
  const size_t order = atoi(option.getValue("-o")[0].c_str());

  // Formulation //
  assemble.start();
  FormulationSteadyWaveScalar<complex<double> > wave(domain, k, order);
  FormulationNeumann neumann(freeSpace, k, order);
  //FormulationNeumann neumann2(outerSpace, k, order);

  // System //
  System<complex<double> > sys(wave);
  SystemHelper<complex<double> >::dirichlet(sys, source, fSourceScal);

  cout << "Free Space (Order: "  << order
       << " --- Wavenumber: "    << k
       << "): " << sys.getSize() << endl;

  // Assemble and Neumann//
  sys.assemble();
  sys.addBorderTerm(neumann);
  //sys.addBorderTerm(neumann2);
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

  // Interpolate //
  if(option.getValue("-interp").size() != 0){
    Mesh           visuMsh(option.getValue("-interp")[0]);
    GroupOfElement visu(visuMsh.getFromPhysical(7));

    fullMatrix<double> point;
    set<const MVertex*, VertexComparator> vertex;
    visu.getAllVertex(vertex);

    {
      set<const MVertex*, VertexComparator>::iterator it  = vertex.begin();
      set<const MVertex*, VertexComparator>::iterator end = vertex.end();

      point.resize(vertex.size(), 3);

      for(size_t i = 0; it != end; it++, i++){
        point(i, 0) = (*it)->x();
        point(i, 1) = (*it)->y();
        point(i, 2) = (*it)->z();
      }
    }

    fullVector<complex<double> > sol;
    sys.getSolution(sol, 0);

    fullMatrix<complex<double> > values;

    Interpolator<complex<double> >::interpolate(sys.getFunctionSpace(),
                                                sys.getDofManager(),
                                                sol, point, values);

    map<const MVertex*, complex<double> > nodeData;

    {
      set<const MVertex*, VertexComparator>::iterator it  = vertex.begin();
      set<const MVertex*, VertexComparator>::iterator end = vertex.end();

      for(size_t i = 0; it != end; it++, i++)
        nodeData.insert
          (pair<const MVertex*, complex<double> >(*it, values(i, 0)));
    }

    {
      NodeSolution nodeSolution;
      nodeSolution.addNodeValue(0, 0, visuMsh, nodeData);
      nodeSolution.write("free_interp");
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
