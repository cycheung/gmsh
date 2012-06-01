#include "Solver.h"

extern "C"{
#include <clapack.h>
}

void Solver::solve(fullMatrix<double>& A, 
		   fullVector<double>& x,
		   fullVector<double>& b){
  
  A.luSolve(b, x);
  /*
  int *ipiv = new int[x.N];
  
  clapack_dgesv(CblasRowMajor, x.N, 1, A.matrix, A.nCol, ipiv, x.v, x.N);

  delete[] ipiv;
  */
}
