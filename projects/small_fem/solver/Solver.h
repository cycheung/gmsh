#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "fullMatrix.h"

/**
   @class Solver
   @brief A linear system solver

   This class handles a solver for the @em linear system
   @c A @c x @c = @c b.

   @note
   The solver is called by the
   @em static method Solver::solve.@n
   So, their is @em no @em need to instantiate a Solver
*/

class Solver{
 public:
  static void solve(fullMatrix<double>& A,
		    fullVector<double>& x,
		    fullVector<double>& b);
};

/**
   @fn Solver::solve
   @param A The Matrix of the system to solve
   @param x The Vector with the futur solution
   ot the system to solve
   @param b The Vector with the
   Right Hand Side of the system to solve

   Solves the linear System Ax = b
   **
 */

#endif
