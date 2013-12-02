#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "fullMatrix.h"
#include "SolverVector.h"
#include "SolverMatrix.h"

/**
   @interface Solver
   @brief Common interface for solvers

   This class is a common interface for linear system solvers.

   A Solver may be of the following scalar types:
   @li Real
   @li Complex

 */

template<typename scalar>
class Solver{
 public:
  virtual ~Solver(void);
  virtual void solve(SolverMatrix<scalar>& A,
                     SolverVector<scalar>& rhs,
                     fullVector<scalar>& x) = 0;
 protected:
  Solver(void);
};

/**
   @fn Solver::~Solver
   Deletes this Solver
   **

   @fn Solver::solve
   @param A A SolverMatrix
   @param rhs A SolverVector
   @param x A fullVector

   Solver the linear system Ax = rhs.

   The matrix and the right hand side can be modified,
   depending on the actual implementation.
   This should be found in the description of the actual solver.
 */


//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SolverInclusion.h"

#endif
