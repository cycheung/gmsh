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
 */

template<typename scalar>
class Formulation{
 public:
  virtual ~Formulation(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId)  const = 0;
  virtual scalar weakB(size_t dofI, size_t dofJ, size_t elementId) const = 0;
  virtual scalar rhs(size_t equationI, size_t elementId)           const = 0;

  virtual bool          isGeneral(void) const = 0;
  virtual const FunctionSpace& fs(void) const = 0;
};

/**
   @fn Formulation::~Formulation
   Deletes this Formulation
   **

   @fn Formulation::weak
   @param dofI The first index of the formulation term
   @param dofJ The second index of the formulation term
   @param elementId The element ID associated with the formulation term
   @return The value of the requested formulation term
   **

   @fn Formulation::weakB
   @param dofI The first index of the formulation term
   @param dofJ The second index of the formulation term
   @param elementId The element ID associated with the formulation term

   This method is only valid when Formulation::isGeneral() is true

   @return The value of the requested second formulation term
   **

   @fn Formulation::rhs
   @param equationI The ith equation of the formulation
   @param elementId The element ID associated
   with the ith equation of the formulation
   @return The value of the ith equation right hand side
   **

   @fn Formulation::isGeneral
   @return Returns
   @li true, if the problem is a generalized eigenvalue problem
   @li false, if not
   **

   @fn Formulation::fs
   @return Returns the FunctionSpace used by this Formulation
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationInclusion.h"

#endif
