#ifndef _FORMULATIONSTEADYWAVE_H_
#define _FORMULATIONSTEADYWAVE_H_

#include <vector>

#include "FunctionSpaceEdge.h"
#include "Formulation.h"

/**
   @class FormulationSteadyWave
   @brief Formulation for the Steady Wave problem

   Formulation for the @em Steady @em Wave problem.

   @todo
   Remove ALL const_cast%S by correcting MElement constness@n
   Allow Hybrid Mesh
 */

class FormulationSteadyWave: public Formulation{
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

  // Function Space //
  FunctionSpaceEdge* fspace;

 public:
  FormulationSteadyWave(const GroupOfElement& goe,
			double k,
			unsigned int order);

  virtual ~FormulationSteadyWave(void);

  virtual double weak(int dofI, int dofJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationSteadyWave::FormulationSteadyWave
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWave of the given 
   @em order and @em wave @em number (@c k)@n

   The given GroupOfElement will be used as the 
   geomtrical @em domain
   **

   @fn FormulationSteadyWave::~FormulationSteadyWave
   Deletes this FormulationSteadyWave
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline double FormulationSteadyWave::rhs(int equationI,
					 const GroupOfDof& god) const{
  return 0;
}

inline const FunctionSpace& FormulationSteadyWave::fs(void) const{
  return *fspace;
}

#endif
