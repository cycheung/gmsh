#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "Matrix.h"
#include "Vector.h"

/**
   @class Solver
   @brief A linear system solver

   This class handles a solver for the @em linear system
   @c A @c x @c = @c b.
 
   @note
   The solver is called by the 
   @em static method Solver::solve.@n
   So, their is @em no @em need to instantiate a Solver

   @todo
   Return the solution in a sperate vector, to save the RHS@n
   Use other stuff than LAPACK
*/

class Solver{
 public:
  static void solve(Matrix& A, Vector<double>& x);
};

/**
   @fn Solver::solve
   @param A The Matrix of the system to solve
   @param x The Vector with the 
   Right Hand Side of the system to solve
   @returns The Vector @c x becomes the solution of the system
 */

#endif
