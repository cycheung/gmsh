#ifndef _FORMULATIONSTEADYWAVESCALAR_H_
#define _FORMULATIONSTEADYWAVESCALAR_H_

#include "FunctionSpaceScalar.h"

#include "TermGradGrad.h"
#include "TermFieldField.h"

#include "Formulation.h"

/**
   @class FormulationSteadyWaveScalar
   @brief Scalar Formulation for the Steady Wave problem

   Scalar Formulation for the steady wave problem
 */

template<typename scalar>
class FormulationSteadyWaveScalar: public Formulation<scalar>{
 private:
  // Speed of medium squared //
  static const double cSquare;

  // Pulsation Squared //
  double omegaSquare;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;

  // Local Terms //
  TermGradGrad*   localTerms1;
  TermFieldField* localTerms2;

 public:
  FormulationSteadyWaveScalar(GroupOfElement& goe,
                              double omega,
                              size_t order);

  virtual ~FormulationSteadyWaveScalar(void);

  virtual bool isGeneral(void) const;

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId)  const;
  virtual scalar weakB(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationSteadyWaveScalar::FormulationSteadyWaveScalar
   @param goe A GroupOfElement
   @param omega A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWaveScalar of the given
   order and pulsation (omega)@n

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationSteadyWaveScalar::~FormulationSteadyWaveScalar
   Deletes this FormulationSteadyWaveScalar
*/

#endif
