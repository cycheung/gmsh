#ifndef _FORMULATION_H_
#define _FORMULATION_H_

#include "GroupOfDof.h"
#include "FunctionSpace.h"

/**
   @interface Formulation
   @brief Base interface of a finite element formulation

   This is the base interface of a finite element formulation.@n

   A Formulation is defined by a @em bilinear @em weak formulation,
   and a @em right @em hand @em side.@n

   @warning
   A formulation is defined @em only on @em GroupOfDof%s.

   @todo
   Add quadrature laws as a paramaeter of a Formulation@n
 */

class Formulation{
 public:
  virtual ~Formulation(void);
  
  virtual double weak(int entityI, int entityJ,
		      const GroupOfDof& god) const = 0;
  
  virtual double rhs(int equationI, 
		     const GroupOfDof& god) const = 0;

  virtual FunctionSpace& fs(void) const = 0;
};

/**
   @fn Formulation::weak
   @param entityI The @em first index of a formulation term 
   @param entityJ The @em second index of the formulation term
   @param god The @em GroupOfDof associated with the formulation term   
   @return The value of the given formulation term

   @fn Formulation::rhs
   @param equationI The @em ith equation of the formulation
   @param god The @em GroupOfDof associated 
   with the @em ith  equation of the formulation
   @return The value of the @em ith equation Right Hand Side
*/

#endif
