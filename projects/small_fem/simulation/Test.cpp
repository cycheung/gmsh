#include <complex>
#include <iostream>

#include "SmallFem.h"

#include "Timer.h"

#include "LineReferenceSpace.h"
#include "TriReferenceSpace.h"
#include "QuadReferenceSpace.h"
#include "TetReferenceSpace.h"
#include "HexReferenceSpace.h"
#include "PyrReferenceSpace.h"
#include "PriReferenceSpace.h"

#include "BasisGenerator.h"
#include "TriLagrangeBasis.h"
#include "LineNodeBasis.h"
#include "LineEdgeBasis.h"
#include "LineNedelecBasis.h"
#include "TriNodeBasis.h"
#include "QuadNedelecBasis.h"

#include "System.h"
#include "SystemHelper.h"

#include "FormulationNeumann.h"
#include "FormulationSteadyWaveScalar.h"
#include "FormulationProjectionScalar.h"
#include "FormulationProjectionVector.h"

#include "Mesh.h"
#include "fullMatrix.h"
#include "GroupOfJacobian.h"

#include "PermutationTree.h"

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "SolverMUMPS.h"

using namespace std;

complex<double> fSource(fullVector<double>& xyz){
  return complex<double>(1, 0);
}

int main(int argc, char** argv){
  // SmallFEM
  SmallFem::Keywords("-msh,-o,-k");
  SmallFem::Initialize(argc, argv);

  // Options //
  const Options& option = SmallFem::getOptions();

  // Get Domains //
  Mesh msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);
  GroupOfElement source = msh.getFromPhysical(5);
  GroupOfElement wall   = msh.getFromPhysical(6);

  // Get Parameters //
  const double puls  = atof(option.getValue("-k")[0].c_str());
  const size_t order = atoi(option.getValue("-o")[0].c_str());

  // Formulations //
  FormulationSteadyWaveScalar<complex<double> > wave(domain, puls, order);
  FormulationNeumann                           neumann(wall, puls, order);

  // System //
  // Init
  System<complex<double> > system(wave);

  // Dirichlet
  SystemHelper<complex<double> >::dirichlet(system, source, fSource);

  // Assemble with Dirichlet
  system.assemble();

  // Assemble Neumann term
  system.addBorderTerm(neumann);

  // Solve
  system.solve();

  // Write Solution //
  FEMSolution<complex<double> > feSol;
  system.getSolution(feSol);
  feSol.write("neumann");

  // Finalize //
  SmallFem::Finalize();
}
