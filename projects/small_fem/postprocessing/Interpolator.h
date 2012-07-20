#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include "GroupOfElement.h"

/**
   @class Interpolator
   @brief Interpolates Entity values on Node%s
   
   This @em class is the mother (by @em inheritence) 
   of all Interpolator%s.@n

   An Interpolator interpolates Entity values 
   (of a GroupOfElement) on the @em Node%s 
   (of the same GroupOfElement).@n

   The Entity @em type to concider is given
   by the @em different @em implementations 
   of this class.
 
   @note
   An Interpolator can be of @em two types
   @li InterpolatorScalar, for @em scalar fiels
   @li InterpolatorVector, for @em vectorial fields
*/

class Interpolator{
 protected:
  bool scalar;
  bool gotInterpolation;

 public:
  virtual ~Interpolator(void);

  virtual void interpolate(const GroupOfElement& mesh) = 0;

  bool isScalar(void) const;

 protected:
  Interpolator(void);
};

/**
   @fn Interpolator::~Interpolator
   @return Deletes the Interpolator

   @fn Interpolator::interpolate
   @param mesh A GroupOfElement to interpolate on
   @return Computes the interpolation of the
   Entity of the given GroupOfElement on the Node%s
   (of the same GroupOfElement)

   @fn Interpolator::isScalar
   @return Returns:
   @li @c true, if this is an Interpolator for 
   @em scalar fields
   @li @c false, if this is an Interpolator for 
   @em vectorial fields
 */ 

//////////////////////
// Inline Functions //
//////////////////////

inline bool Interpolator::isScalar(void) const{
  return scalar;
}

#endif
