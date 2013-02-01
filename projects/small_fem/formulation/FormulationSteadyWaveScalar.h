#ifndef _FORMULATIONSTEADYWAVESCALAR_H_
#define _FORMULATIONSTEADYWAVESCALAR_H_

#include "FunctionSpaceScalar.h"

#include "TermGradGrad.h"
#include "TermFieldField.h"

#include "Formulation.h"

/**
   @class FormulationSteadyWaveScalar
   @brief Scalar Formulation for the Steady Wave problem

   Scalar Formulation for the @em Steady @em Wave problem
 */

class FormulationSteadyWaveScalar: public Formulation{
 private:
  // Physical Values //
  static const double mu;
  static const double eps;

  // Wave Number Squared //
  double kSquare;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;

  // Local Terms //
  TermGradGrad*   localTerms1;
  TermFieldField* localTerms2;

 public:
  FormulationSteadyWaveScalar(GroupOfElement& goe,
			      double k,
			      unsigned int order);

  virtual ~FormulationSteadyWaveScalar(void);

  virtual double weak(unsigned int dofI, unsigned int dofJ,
		      const GroupOfDof& god) const;

  virtual double rhs(unsigned int equationI,
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

#endif
