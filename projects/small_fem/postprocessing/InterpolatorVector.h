#ifndef _INTERPOLATORVECTOR_H_
#define _INTERPOLATORVECTOR_H_

#include <vector>
#include "Vector.h"
#include "Interpolator.h"

/**
   @class InterpolatorVector
   @brief Interpolator for vectorial fields
   
   This an Interpolator for @em vectorial fields.

   @note
   An InterpolatorVector can't be instantiate (but descendant can)
*/

class InterpolatorVector: public Interpolator{
 protected:
  std::vector<Vector<double>*>* nodeValue;

 public:
  virtual ~InterpolatorVector(void);

  std::vector<Vector<double>*>* getNodeValue(void) const;

 protected:
  InterpolatorVector(void);
};

/**
   @fn InterpolatorVector::~InterpolatorVector
   @return Deletes the InterpolatorVector

   @fn InterpolatorVector::getNodeValue
   @return Returns the interpolated 
   @em vectorial field
 */


#endif
