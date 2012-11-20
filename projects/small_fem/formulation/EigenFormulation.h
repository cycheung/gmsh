#ifndef _EIGENFORMULATION_H_
#define _EIGENFORMULATION_H_

#include "GroupOfDof.h"
#include "FunctionSpace.h"

/**
   @interface EigenFormulation
   @brief Base interface of a Finite Element 
   Eigenvalue Problem formulation

   This is the base interface of a Finite Element 
   Eigenvalue formulation.@n

   A Finite Element Eigenvalue Problem can be of @em two types:
   @li Generalized EigenFormulation, if the problem 
   is of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{B})\mathbf{x} = \mathbf{b}@f$
   @li Not Generalized EigenFormulation, if the problem 
   is of the type
   @f$\qquad(\mathbf{A} - \lambda{}\mathbf{I})\mathbf{x} = \mathbf{b}@f$

   A Eigenvalue formulation is defined by:
   @li A Function Space 
   @see FunctionSpace class
   @li Two @em or one Bilinear Weak formulation 
   (depending if the Problem is a @em generalized Eigenvalue problem or not).@n

   @warning
   A EigenFormulation is defined @em only on @em GroupOfDof%s.

   @todo
   Add quadrature laws as a paramaeter of a EigenFormulation@n
   Allow evaluation of non GroupOfDof related Dof.
 */

class EigenFormulation{
 public:
  virtual ~EigenFormulation(void);
  
  virtual bool isGeneral(void) const = 0;

  virtual double weakA(int dofI, int dofJ,
		       const GroupOfDof& god) const = 0;

  virtual double weakB(int dofI, int dofJ,
		       const GroupOfDof& god) const = 0;
  
  virtual const FunctionSpace& fs(void) const = 0;
};

/**
   @fn EigenFormulation::~EigenFormulation
   Deletes this EigenFormualtion
   **

   @fn EigenFormulation::isGeneral
   @return Returns 
   @li @c true, if the problem is a @em generalized Eigenvalue problem
   @li @c false, if not
   **

   @fn EigenFormulation::weakA
   @param dofI The @em first index of a Eigenvalue formulation term 
   @param dofJ The @em second index of the Eigenvalue formulation term
   @param god The @em GroupOfDof associated with the Eigenvalue formulation term   
   @return The value of the requested Eigenvalue formulation term, 
   @em for the @f$A@f$ Matrix
   **

   @fn EigenFormulation::weakB
   @param dofI The @em first index of a Eigenvalue formulation term 
   @param dofJ The @em second index of the Eigenvalue formulation term
   @param god The @em GroupOfDof associated with the Eigenvalue formulation term   
   @return The value of the requested Eigenvalue formulation term, 
   @em for the @f$B@f$ Matrix
   
   @warning This method has a meaning only if the Matrix @f$B@f$ exists@n
   That is if EigenFormulation::isGeneral() is @c true
   **

   @fn EigenFormulation::fs
   @return Returns the Function Space used by this EigenFormulation
   **
*/

#endif
