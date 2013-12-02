#ifndef _FORMULATION_H_
#define _FORMULATION_H_

#include <string>
#include "FunctionSpace.h"

/**
   @interface Formulation
   @brief Base interface of a finite element formulation

   This is the base interface of a finite element formulation.

   A finite element problem is of the type:
   @f$\mathbf{A}~\mathbf{x} = \mathbf{b}@f$

   A finite element eigenvalue problem is a of the type:
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$

   A generalized finite element eigenvalue problem is a of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$

   To conclude a Formulation may be of multiple scalar type.
   See FormulationTyped for more informations.
 */

class Formulation{
 public:
  virtual ~Formulation(void);

  virtual bool                 isGeneral(void) const = 0;
  virtual const FunctionSpace& fs(void)        const = 0;
  virtual std::string          getType(void)   const = 0;
};

/**
   @fn Formulation::~Formulation
   Deletes this Formulation
   **

   @fn Formulation::isGeneral
   @return Returns
   @li true, if the problem is a generalized eigenvalue problem
   @li false, if not
   **

   @fn Formulation::fs
   @return Returns the FunctionSpace used by this Formulation
   **

   @fn Formulation::getType
   @return Returns a string specifying of which type this Formulation is
   @see FormulationTyped
*/

#endif
