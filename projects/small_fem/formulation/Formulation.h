#ifndef _FORMULATION_H_
#define _FORMULATION_H_

#include "GroupOfDof.h"
#include "FunctionSpace.h"

/**
   @interface Formulation
   @brief Base interface of a Finite Element formulation

   This is the base interface of a Finite Element formulation.@n

   A Finite Element Problem is of the type:@n
   @f$\mathbf{A}~\mathbf{x} = \mathbf{b}@f$@n

   A Formulation is defined by:
   @li A Function Space
   @see FunctionSpace class
   @li A Bilinear Weak Formulation
   @li A Right hand Side.@n

   @warning
   A formulation is defined @em only on @em GroupOfDof%s.

   @todo
   Add quadrature laws as a paramaeter of a Formulation@n
   Allow evaluation of non GroupOfDof related Dof.@n
   virtual weak() is not a good idea with 'fast' assembly
 */

class Formulation{
 public:
  virtual ~Formulation(void);

  virtual bool isGeneral(void) const = 0;

  virtual double weak(unsigned int dofI, unsigned int dofJ,
		      const GroupOfDof& god) const = 0;

  virtual double weakB(unsigned int dofI, unsigned int dofJ,
		       const GroupOfDof& god) const = 0;

  virtual double rhs(unsigned int equationI,
		     const GroupOfDof& god) const = 0;

  virtual const FunctionSpace& fs(void) const = 0;
};

/**
   @fn Formulation::~Formulation
   Deletes this Formulation
   **

   @fn Formulation::isGeneral
   @return Returns
   @li @c true, if the problem is a @em generalized Eigenvalue problem
   @li @c false, if not
   **

   @fn Formulation::weak
   @param dofI The @em first index of the formulation term
   @param dofJ The @em second index of the formulation term
   @param god The @em GroupOfDof associated with the formulation term
   @return The value of the requested formulation term
   **

   @fn Formulation::weakB
   @param dofI The @em first index of the formulation term
   @param dofJ The @em second index of the formulation term
   @param god The @em GroupOfDof associated with the formulation term

   This method is only valid when Formulation::isGeneral() is @c true

   @return The value of the requested @em second formulation term
   **

   @fn Formulation::rhs
   @param equationI The @em ith equation of the formulation
   @param god The @em GroupOfDof associated
   with the @em ith  equation of the formulation
   @return The value of the @em ith equation Right Hand Side
   **

   @fn Formulation::fs
   @return Returns the Function Space used by this Formulation
*/

#endif
