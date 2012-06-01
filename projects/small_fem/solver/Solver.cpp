#include "Solver.h"

extern "C"{
#include <clapack.h>
}

void Solver::solve(Matrix& A, Vector<double>& x){
  int *ipiv = new int[x.N];
  
  clapack_dgesv(CblasRowMajor, x.N, 1, A.matrix, A.nCol, ipiv, x.v, x.N);

  delete[] ipiv;
}
