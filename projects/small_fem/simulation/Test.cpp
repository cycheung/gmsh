#include <iostream>

#include "SmallFem.h"

#include "Timer.h"
#include "TriLagrangeBasis.h"
#include "LineReferenceSpace.h"
#include "TriReferenceSpace.h"
#include "QuadReferenceSpace.h"
#include "TetReferenceSpace.h"
#include "HexReferenceSpace.h"
#include "PyrReferenceSpace.h"
#include "PriReferenceSpace.h"

#include "LineNodeBasis.h"
#include "LineEdgeBasis.h"
#include "LineNedelecBasis.h"
#include "TriNodeBasis.h"
#include "QuadNedelecBasis.h"

#include "Mesh.h"
#include "fullMatrix.h"
#include "GroupOfJacobian.h"

#include "PermutationTree.h"

#include "SolverMatrix.h"
#include "SolverVector.h"
#include "SolverMUMPS.h"

using namespace std;

int main(int argc, char** argv){
  SmallFem::Initialize(argc, argv);

  SolverMatrix<complex<double> > a(5, 5);
  SolverVector<complex<double> > b(5);
  fullVector<complex<double> >   x;
  SolverMUMPS<complex<double> >  solver;

  a.add(1 - 1, 2 - 1, complex<double>(+3.0, 1.0));
  a.add(2 - 1, 3 - 1, complex<double>(-3.0, 1.0));
  a.add(4 - 1, 3 - 1, complex<double>(+2.0, 1.0));
  a.add(5 - 1, 5 - 1, complex<double>(+1.0, 1.0));
  a.add(2 - 1, 1 - 1, complex<double>(+3.0, 1.0));
  a.add(1 - 1, 1 - 1, complex<double>(+2.0, 1.0));
  a.add(5 - 1, 2 - 1, complex<double>(+4.0, 1.0));
  a.add(3 - 1, 4 - 1, complex<double>(+2.0, 1.0));
  a.add(2 - 1, 5 - 1, complex<double>(+6.0, 1.0));
  a.add(3 - 1, 2 - 1, complex<double>(-1.0, 1.0));
  a.add(1 - 1, 3 - 1, complex<double>(+4.0, 1.0));
  a.add(3 - 1, 3 - 1, complex<double>(+1.0, 1.0));

  b.add(1 - 1, complex<double>(+20.0, 0.0));
  b.add(2 - 1, complex<double>(+24.0, 0.0));
  b.add(3 - 1, complex<double>(+09.0, 0.0));
  b.add(4 - 1, complex<double>(+06.0, 0.0));
  b.add(5 - 1, complex<double>(+13.0, 0.0));

  solver.solve(a, b, x);

  for(int i = 0; i < x.size(); i++)
    cout << x(i) << endl;

  a.writeToMatlabFile("a.m", "a");

  SmallFem::Finalize();
}
