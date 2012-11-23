#ifndef _FORMULATIONSTEADYWAVEVECTOR_H_
#define _FORMULATIONSTEADYWAVEVECTOR_H_

#include <vector>

#include "FunctionSpaceEdge.h"
#include "Formulation.h"

/**
   @class FormulationSteadyWaveVector
   @brief Vectorial Formulation for the Steady Wave problem

   Vectorial Formulation for the @em Steady @em Wave problem.

   @todo
   Remove ALL const_cast%S by correcting MElement constness@n
   Allow Hybrid Mesh
 */

class FormulationSteadyWaveVector: public Formulation{
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
  FormulationSteadyWaveVector(const GroupOfElement& goe,
			      double k,
			      unsigned int order);

  virtual ~FormulationSteadyWaveVector(void);

  virtual double weak(int dofI, int dofJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationSteadyWaveVector::FormulationSteadyWaveVector
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWaveVector of the given 
   @em order and @em wave @em number (@c k)@n

   The given GroupOfElement will be used as the 
   geomtrical @em domain
   **

   @fn FormulationSteadyWaveVector::~FormulationSteadyWaveVector
   Deletes this FormulationSteadyWaveVector
*/

//////////////////////
// Inline Functions //
//////////////////////

inline double FormulationSteadyWaveVector::rhs(int equationI,
					       const GroupOfDof& god) const{
  return 0;
}

inline const FunctionSpace& FormulationSteadyWaveVector::fs(void) const{
  return *fspace;
}

#endif
