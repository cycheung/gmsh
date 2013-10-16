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

#include "TriNodeBasis.h"
#include "LineNodeBasis.h"

#include "Mesh.h"
#include "fullMatrix.h"
#include "GroupOfJacobian.h"

#include "PermutationTree.h"

#include "SparseMatrix.h"
#include "SolverMUMPS.h"

using namespace std;

int main(int argc, char** argv){
  //SmallFem::Initialize(argc, argv);

  size_t n = 2;

  SparseMatrix A(n, n);
  fullVector<double> b(n);
  fullVector<double> x(n);

  A.add(0, 0, 1);
  A.add(1, 1, 2);

  b(0) = 1;
  b(1) = 4;

  SolverMUMPS solver;
  solver.solve(A, b, x);

  for(size_t i = 0; i < n; i++)
    cout << x(i) << endl;

  //SmallFem::Finalize();

  /*
  Timer time;
  time.start();

  PriReferenceSpace ref;
  cout << ref.toString() << endl;

  time.stop();

  cout << "Time: " << time.time() << " " << time.unit() << endl;
  */
  /*
  SmallFem::Initialize(argc, argv);
  #pragma omp parallel for
  for(int i = 0; i < 9; i++)
    printf("%d\n", i);

  QuadReferenceSpace ref;
  cout << ref.toString() << endl;

  SmallFem::Finalize();
  */
  return 0;

  /*
  const size_t N = 387 * 10;
  const size_t M = 28224 * 10;

  Timer   timer;
  double* ma = new double[N * M];

  for(size_t i = 0; i < N * M; i++)
    ma[i] = 0;

  timer.start();

#pragma omp parallel for
  for(size_t j = 0; j < M; j++)
    for(size_t i = 0; i < N; i++)
      ma[i + N * j] = i + N * j;

  timer.stop();

  cout << timer.time() << " "
       << timer.unit() << endl
       << flush;

  delete[] ma;
  return 0;
  */
}
