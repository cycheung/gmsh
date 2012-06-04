#ifndef _INTERPOLATORSCALAR_H_
#define _INTERPOLATORSCALAR_H_

#include <vector>
#include "Interpolator.h"

/**
   @class InterpolatorScalar
   @brief Interpolator for scalar fields
   
   This an Interpolator for @em scalar fields.
   
   @note
   An InterpolatorScalar can't be instantiate (but descendant can)
*/

class InterpolatorScalar: public Interpolator{
 protected:
  std::vector<double>* nodeValue;

 public:
  virtual ~InterpolatorScalar(void);

  const std::vector<double>& getNodeValue(void) const;

 protected:
  InterpolatorScalar(void);
};

/**
   @fn InterpolatorScalar::~InterpolatorScalar
   @return Deletes the InterpolatorScalar

   @fn InterpolatorScalar::getNodeValue
   @return Returns the interpolated 
   @em scalar field
 */

#endif
