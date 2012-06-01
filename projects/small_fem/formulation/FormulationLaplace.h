#ifndef _FORMULATIONLAPLACE_H_
#define _FORMULATIONLAPLACE_H_

#include "Formulation.h"
#include "Polynomial.h"
#include "Vector.h"
#include "TriNodeBasis.h"
#include "InterpolatorNode.h"

/**
   @class FormulationLaplace
   @brief Formulation for the Laplace problem

   Formulation for the @em Laplace problem.
 */

class FormulationLaplace: public Formulation{
 private:
  // Gaussian Quadrature Data //
  int G;
  double gx[4];
  double gy[4];
  double gw[4];

  // Basis //
  TriNodeBasis*       base;
  Vector<Polynomial>* gradBasis;
  int                 basisSize;

  // Interpolator //
  InterpolatorNode* interp;

 public:
  FormulationLaplace(void);

  virtual ~FormulationLaplace(void);

  virtual double weak(const int nodeI, const int nodeJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(const int equationI,
		     const GroupOfDof& god) const;

  virtual Interpolator& interpolator(void) const;
};

/**
   @fn FormulationLaplace::FormulationLaplace
   @return Returns a new FormulationLaplace
 
   @fn FormulationLaplace::~FormulationLaplace
   @return Deletes this FormulationLaplace
*/

//////////////////////
// Inline Functions //
//////////////////////

inline double FormulationLaplace::rhs(const int equationI,
				      const GroupOfDof& god) const{
  return 0;
}

inline Interpolator& FormulationLaplace::interpolator(void) const{
  return *interp;
}

#endif
