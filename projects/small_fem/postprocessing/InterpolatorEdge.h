#ifndef _INTERPOLATOREDGE_H_
#define _INTERPOLATOREDGE_H_

#include <vector>
#include "InterpolatorVector.h"
#include "BasisVector.h"
#include "Polynomial.h"

/**
   @class InterpolatorEdge
   @brief Interpolator for Edge%s

   This class is a @em vectorial Interpolator
   for values defined on the @em Edge%s of a 
   given Mesh.
 */

class InterpolatorEdge: public InterpolatorVector{
 private:
  const std::vector<Polynomial>* basis;
  int bSize;

  const Mesh* msh;  
  int nNode;

  std::vector<bool>* isInterpolated;

 public:
  InterpolatorEdge(const BasisVector& basis);
   
  virtual ~InterpolatorEdge(void);
  
  virtual void interpolate(const Mesh& mesh);
  
 private:
  void interpolateEdgeElement(void);
};

/**
   @fn InterpolatorEdge::InterpolatorEdge
   @param basis The Basis to use for 
   interpolation
   @return Returns a new InterpolatorEdge

   @fn InterpolatorEdge::~InterpolatorEdge
   @return Deletes this InterpolatorEdge
 */

#endif
