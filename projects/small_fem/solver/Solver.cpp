#include "Solver.h"

void Solver::solve(fullMatrix<double>& A,
		   fullVector<double>& x,
		   fullVector<double>& b){

  A.luSolve(b, x);
}
