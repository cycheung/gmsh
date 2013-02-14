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
		       unsigned int order);

  virtual ~FormulationVibration(void);

  virtual bool isGeneral(void) const;

  virtual double weak(unsigned int dofI, unsigned int dofJ,
                      const GroupOfDof& god) const;

  virtual double weakB(unsigned int dofI, unsigned int dofJ,
		       const GroupOfDof& god) const;

  virtual double rhs(unsigned int equationI,
		     const GroupOfDof& god) const;

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
