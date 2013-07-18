#ifndef _FORMULATIONVIBRATION_H_
#define _FORMULATIONVIBRATION_H_

#include "FunctionSpaceScalar.h"
#include "TermGradGrad.h"
#include "Formulation.h"

/**
   @class FormulationVibration
   @brief Formulation for the Vibration problem

   Formulation for the @em Vibration problem
 */

class FormulationVibration: public Formulation{
 private:
  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;

  // Local Terms //
  TermGradGrad* localTerms;

 public:
  FormulationVibration(GroupOfElement& goe,
                       size_t order);

  virtual ~FormulationVibration(void);

  virtual bool isGeneral(void) const;

  virtual double weak(size_t dofI, size_t dofJ,
                      size_t elementId) const;

  virtual double weakB(size_t dofI, size_t dofJ,
                       size_t elementId) const;

  virtual double rhs(size_t equationI,
                     size_t elementId) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationVibration::FormulationVibration
   @param goe A GroupOfElement
   @param order A natural number

   Instantiates a new FormulationVibration of the given order@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationVibration::~FormulationVibration
   Deletes this FormulationVibration
*/

#endif
