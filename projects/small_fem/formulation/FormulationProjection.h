#ifndef _FORMULATIONPROJECTION_H_
#define _FORMULATIONPROJECTION_H_

#include "Formulation.h"
#include "Vector.h"
#include "Polynomial.h"
#include "TriNedelecBasis.h"
#include "InterpolatorEdge.h"

/**
   @class FormulationProjection
   @brief Formulation for the Projection problem

   Vectorial Formulation for the @em Projection problem.
 */

class FormulationProjection: public Formulation{
 private:
  // Gaussian Quadrature Data //
  int G;
  double gx[4];
  double gy[4];
  double gw[4];

  // Basis //
  TriNedelecBasis*          baseGen;
  const Vector<Polynomial>* basis;
  int                       basisSize;

  // Vector to Project //
  Vector<double>* f;

  // Interpolator //
  InterpolatorEdge* interp;

 public:
  FormulationProjection(Vector<double>& vectorToProject);
  
  virtual ~FormulationProjection(void);

  virtual double weak(const int edgeI, const int edgeJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(const int equationI,
		     const GroupOfDof& god) const;

  virtual Interpolator& interpolator(void) const;
};

/**
   @fn FormulationProjection::FormulationProjection
   @param vectorToProject A Vector<double>
   @return Returns a new FormulationProjection to project
   the given Vector
 
   @fn FormulationProjection::~FormulationProjection
   @return Deletes the this FormulationProjection
*/

//////////////////////
// Inline Functions //
//////////////////////

inline Interpolator& FormulationProjection::interpolator(void) const{
  return *interp;
}


#endif
