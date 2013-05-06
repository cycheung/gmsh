#include <iostream>
#include "Timer.h"

using namespace std;

int main(int argc, char** argv){
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
}
