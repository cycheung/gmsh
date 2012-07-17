#ifndef _FORMULATION_H_
#define _FORMULATION_H_

#include "GroupOfDof.h"
#include "Interpolator.h"
#include "Dof.h"

/**
   @interface Formulation
   @brief Base interface of a finite element formulation

   This is the base interface of a finite element formulation.@n

   A Formulation is defined by a @em bilinear @em weak formulation,
   and a @em right @em hand @em side.@n

   A Formulation shall also give an Interpolator, for
   handling solutions of this Formulation.
   
   @warning
   A formulation is defined @em only on @em GroupOfDof%s.

   @todo
   Add quadrature laws as a paramaeter of a Formulation@n
   Remove dependance on GroupOfDof%s
 */

class Formulation{
 public:
  virtual ~Formulation(void);
  
  virtual double weak(const int entityI, const int entityJ,
		      const GroupOfDof& god) const = 0;
  
  virtual double rhs(const int equationI,
		     const GroupOfDof& god) const = 0;

  virtual const std::vector<Dof*> getAllDofs(void) const = 0;
  
  //virtual Interpolator& interpolator(void) const = 0;
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

   @fn Formulation::getAllDofs
   @return Returns all the Dofs that are part of 
   the formulation
   @note Note that the first call to this method
   will also generate the requested Dofs

   @fn Formulation::interpolator
   @return Returns the Interpolator associated with
   this Formulation
*/

//////////////////////
// Inline Functions //
//////////////////////

inline Formulation::~Formulation(void){
}

#endif
