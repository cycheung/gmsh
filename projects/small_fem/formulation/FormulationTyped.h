#ifndef _FORMULATIONTYPED_H_
#define _FORMULATIONTYPED_H_

#include "Formulation.h"

/**
   @interface FormultationTyped
   @brief Type of scalar used by a finite element formulation

   This interface adds to Forumlation the notion of scalar type.
   A scalar may be:
   @li real
   @li complex

   This interface enables access to the finite element weak formulation.
 */

template<typename scalar>
class FormulationTyped: public Formulation{
 public:
  virtual ~FormulationTyped(void) {};

  virtual scalar weak(size_t dofI, size_t dofJ,
                      size_t elementId) const = 0;

  virtual scalar weakB(size_t dofI, size_t dofJ,
                       size_t elementId) const = 0;

  virtual scalar rhs(size_t equationI,
                     size_t elementId) const = 0;

  virtual std::string getType(void) const;
};


/**
   @fn FormulationTyped::~FormulationTyped
   Deletes this FormulationTyped
   **

   @fn FormulationTyped::weak
   @param dofI The first index of the formulation term
   @param dofJ The second index of the formulation term
   @param elementId The element ID associated with the formulation term
   @return The value of the requested formulation term
   **

   @fn FormulationTyped::weakB
   @param dofI The first index of the formulation term
   @param dofJ The second index of the formulation term
   @param elementId The element ID associated with the formulation term

   This method is only valid when FormulationTyped::isGeneral() is true

   @return The value of the requested second formulation term
   **

   @fn FormulationTyped::rhs
   @param equationI The ith equation of the formulation
   @param elementId The element ID associated
   with the ith equation of the formulation
   @return The value of the ith equation right hand side
 */

#endif
