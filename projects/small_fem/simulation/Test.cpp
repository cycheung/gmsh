#include <iostream>

#include "SmallFem.h"

#include "Timer.h"
#include "TriLagrangeBasis.h"
#include "LineReferenceSpace.h"
#include "TriReferenceSpace.h"
#include "QuadReferenceSpace.h"
#include "TetReferenceSpace.h"

#include "Mesh.h"
#include "fullMatrix.h"
#include "GroupOfJacobian.h"

using namespace std;

int main(int argc, char** argv){
  SmallFem::Initialize(argc, argv);
  /*
  Timer timer;
  timer.start();


  Mesh msh(argv[1]);


  GroupOfElement domain = msh.getFromPhysical(7);

  fullMatrix<double> p(1, 3);
  p(0, 0) = 0;
  p(0, 1) = 0;
  p(0, 2) = 0;

  Timer timerGoj;
  timerGoj.start();
  GroupOfJacobian goj(domain, p, "invert");
  timerGoj.stop();
  timer.stop();

  cout << timerGoj.time() << " " << timerGoj.unit() << endl;
  cout << timer.time() << " " << timer.unit() << endl;
  cout << goj.toString() << endl;
  */
  SmallFem::Finalize();

  return 0;

  /*
  //TriLagrangeBasis basis(atoi(argv[1]));

  TetReferenceSpace ref;
  //cout << ref.toString() << endl;

  return 0;
  */
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
