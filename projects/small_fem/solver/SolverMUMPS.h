#ifndef _SOLVERMUMPS_H_
#define _SOLVERMUMPS_H_

#include <vector>
#include <complex>
#include "Solver.h"
#include "mumps_c_types.h"

/**
   @class SolverMUMPS
   @brief A Solver using MUMPS library

   This class implements a Solver using the
   MUltifrontal Massively Parallel sparse direct %Solver (MUMPS) library.

   This library can be download at
    <a href="http://mumps.enseeiht.fr">http://mumps.enseeiht.fr</a> or
    <a href="http://graal.ens-lyon.fr/MUMPS">http://graal.ens-lyon.fr/MUMPS</a>.
*/

template<typename scalar>
class SolverMUMPS: public Solver<scalar>{
 public:
  SolverMUMPS(void);

  virtual ~SolverMUMPS(void);

  virtual void solve(SolverMatrix<scalar>& A,
                     SolverVector<scalar>& rhs,
                     fullVector<scalar>& x);
 private:
  void toMUMPSComplex(std::vector<std::complex<double> >& data,
                      mumps_double_complex** out);
  void toMUMPSComplex(SolverVector<std::complex<double> >& data,
                      mumps_double_complex** out);
  void fromMUMPSComplex(mumps_double_complex** in,
                        size_t size,
                        std::vector<std::complex<double> >& data);
  void fromMUMPSComplex(mumps_double_complex** in,
                        size_t size,
                        SolverVector<std::complex<double> >& data);
};

/**
   @fn SolverMUMPS::SolverMUMPS
   Instanciates a new SolverMUMPS
   **

   @fn SolverMUMPS::~SolverMUMPS
   Deletes this SolverMUMPS
*/


//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "SolverMUMPSInclusion.h"

#endif
