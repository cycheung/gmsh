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

complex<double> fScal(fullVector<double>& xyz){
  double tmp =
    sin(10 * xyz(0)) +
    sin(10 * xyz(1)) +
    sin(10 * xyz(2));

  return complex<double>(1, 1) * tmp;
}

fullVector<complex<double> > fVect(fullVector<double>& xyz){
  complex<double> tmp = complex<double>(1, 1);
  fullVector<complex<double> > res(3);

  res(0) = sin(10 * xyz(0)) * tmp;
  res(1) = sin(10 * xyz(1)) * tmp;
  res(2) = sin(10 * xyz(2)) * tmp;

  return res;
}

int main(int argc, char** argv){
  // SmallFEM
  SmallFem::Keywords("-msh,-o");
  SmallFem::Initialize(argc, argv);

  // Options //
  const Options& option = SmallFem::getOptions();

  // Input //
  Mesh           msh(option.getValue("-msh")[0]);
  GroupOfElement domain = msh.getFromPhysical(7);
  int            order  = atoi(option.getValue("-o")[0].c_str());

  // Formulation //
  /*
  Basis* basis = BasisGenerator::generate(domain.get(0).getType(),
                                          0, order, "hierarchical");

  FunctionSpaceScalar fSpace(domain, *basis);
  FormulationProjectionScalar<complex<double> > projection(fScal, fSpace);
  */

  Basis* basis = BasisGenerator::generate(domain.get(0).getType(),
                                          1, order, "hierarchical");

  FunctionSpaceVector fSpace(domain, *basis);
  FormulationProjectionVector<complex<double> > projection(fVect, fSpace);

  // System //
  System<complex<double> > sys(projection);
  sys.assemble();
  sys.solve();

  // Solution //
  FEMSolution<complex<double> > feSol;
  sys.getSolution(feSol);
  feSol.write("test");

  // Finalize //
  delete basis;
  SmallFem::Finalize();
}
