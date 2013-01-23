#ifndef _FORMULATIONSTEADYWAVESCALAR_H_
#define _FORMULATIONSTEADYWAVESCALAR_H_

#include <vector>

#include "FunctionSpaceScalar.h"
#include "Formulation.h"

/**
   @class FormulationSteadyWaveScalar
   @brief Scalar Formulation for the Steady Wave problem

   Scalar Formulation for the @em Steady @em Wave problem.

   @todo
   Remove ALL const_cast%S by correcting MElement constness@n
   Allow Hybrid Mesh
 */

class FormulationSteadyWaveScalar: public Formulation{
 private:
  // Physical Values //
  static const double mu;
  static const double eps;

  // Wave Number Squared //
  double kSquare;

  // Gaussian Quadrature Data (Term One) //
  int G1;
  fullMatrix<double>* gC1;
  fullVector<double>* gW1;

  // Gaussian Quadrature Data (Term Two) //
  int G2;
  fullMatrix<double>* gC2;
  fullVector<double>* gW2;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;

 public:
  FormulationSteadyWaveScalar(GroupOfElement& goe,
			      double k,
			      unsigned int order);

  virtual ~FormulationSteadyWaveScalar(void);

  virtual double weak(int dofI, int dofJ,
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationSteadyWaveScalar::FormulationSteadyWaveScalar
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWaveScalar of the given
   @em order and @em wave @em number (@c k)@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationSteadyWaveScalar::~FormulationSteadyWaveScalar
   Deletes this FormulationSteadyWaveScalar
*/

//////////////////////
// Inline Functions //
//////////////////////

inline double FormulationSteadyWaveScalar::rhs(int equationI,
					       const GroupOfDof& god) const{
  return 0;
}

inline const FunctionSpace& FormulationSteadyWaveScalar::fs(void) const{
  return *fspace;
}

#endif
