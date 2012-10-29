#ifndef _EIGENFORMULATION_H_
#define _EIGENFORMULATION_H_

#include "GroupOfDof.h"
#include "FunctionSpace.h"

/**
   @interface EigenFormulation
   @brief Base interface of a Finite Element 
   Eigenvalue Problem Formulation

   This is the base interface of a Finite Element 
   Eigenvalue Formulation.@n

   A Eigenvalue Formulation is defined by:
   @li A Function Space 
   @see FunctionSpace class
   @li @em Two Bilinear Weak Formulation
   @li A Right hand Side.@n

   @warning
   A EigenFormulation is defined @em only on @em GroupOfDof%s.

   @todo
   Add quadrature laws as a paramaeter of a EigenFormulation@n
   Allow evaluation of non GroupOfDof related Dof.
 */

class EigenFormulation{
 public:
  virtual ~EigenFormulation(void);
  
  virtual double weakA(int dofI, int dofJ,
		       const GroupOfDof& god) const = 0;

  virtual double weakB(int dofI, int dofJ,
		       const GroupOfDof& god) const = 0;
  
  virtual double rhs(int equationI, 
		     const GroupOfDof& god) const = 0;

  virtual const FunctionSpace& fs(void) const = 0;
};

/**
   @fn EigenFormulation::~EigenFormulation
   Deletes this EigenFormualtion
   **

   @fn EigenFormulation::weakA
   @param dofI The @em first index of a Eigenvalue Formulation term 
   @param dofJ The @em second index of the Eigenvalue Formulation term
   @param god The @em GroupOfDof associated with the eigenFormulation term   
   @return The value of the given Eigenvalue Formulation term
   **

   @fn EigenFormulation::rhs
   @param equationI The @em ith equation of the eigenFormulation
   @param god The @em GroupOfDof associated 
   with the @em ith  equation of the eigenFormulation
   @return The value of the @em ith equation Right Hand Side
   **

  @fn EigenFormulation::fs
  @return Returns the Function Space used by this EigenFormulation
  **
*/

#endif
