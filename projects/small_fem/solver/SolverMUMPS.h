#ifndef _SOLVERMUMPS_H_
#define _SOLVERMUMPS_H_

#include "Solver.h"

/**
   @class SolverMUMPS
   @brief A Solver using MUMPS library

   This class implements a Solver using the
   MUltifrontal Massively Parallel sparse direct %Solver (MUMPS) library.

   This library can be download at
    <a href="http://mumps.enseeiht.fr">http://mumps.enseeiht.fr</a> or
    <a href="http://graal.ens-lyon.fr/MUMPS">http://graal.ens-lyon.fr/MUMPS</a>.
*/

class SolverMUMPS: public Solver{
 public:
  SolverMUMPS(void);

  virtual ~SolverMUMPS(void);

  virtual void solve(SolverMatrix& A,
                     SolverVector& rhs,
                     fullVector<double>& x);
};

/**
   @fn SolverMUMPS::SolverMUMPS
   Instanciates a new SolverMUMPS
   **

   @fn SolverMUMPS::~SolverMUMPS
   Deletes this SolverMUMPS
*/

#endif
