#ifndef _INTERPOLATORNODE_H_
#define _INTERPOLATORNODE_H_

#include "InterpolatorScalar.h"

/**
   @class InterpolatorNode
   @brief Interpolator for Edge%s

   This class is a @em scalar Interpolator
   for values defined on the @em Node%s of a 
   given Mesh.
 */

class InterpolatorNode: public InterpolatorScalar{
 public:
  InterpolatorNode(void);
   
  virtual ~InterpolatorNode(void);
  
  virtual void interpolate(const Mesh& mesh);
};

/**
   @fn InterpolatorNode::InterpolatorNode
   @return Returns a new InterpolatorNode

   @fn InterpolatorNode::~InterpolatorNode
   @return Deletes this InterpolatorNode
 */

#endif
