#ifndef _FORMULATIONLAPLACE_H_
#define _FORMULATIONLAPLACE_H_

#include <vector>
#include "Formulation.h"
#include "Polynomial.h"
#include "TriNodeBasis.h"
//#include "InterpolatorNode.h"

/**
   @class FormulationLaplace
   @brief Formulation for the Laplace problem

   Formulation for the @em Laplace problem.
 */

class FormulationLaplace: public Formulation{
 private:
  // Gaussian Quadrature Data //
  int G;
  fullMatrix<double>* gC;
  fullVector<double>* gW;

  // Basis //
  TriNodeBasis*            base;
  std::vector<Polynomial>* gradBasis;
  int                      basisSize;

  // Interpolator //
  //InterpolatorNode* interp;

 public:
  FormulationLaplace(void);

  virtual ~FormulationLaplace(void);

  virtual double weak(const int nodeI, const int nodeJ, 
		      const GeoDof& god) const;

  virtual double rhs(const int equationI,
		     const GeoDof& god) const;

  virtual const std::vector<Dof*> getAllDofs(void) const;

  //virtual Interpolator& interpolator(void) const;
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
				      const GeoDof& god) const{
  return 0;
}
/*
inline Interpolator& FormulationLaplace::interpolator(void) const{
  return *interp;
}
*/
inline const std::vector<Dof*> 
FormulationLaplace::getAllDofs(void) const{
  return std::vector<Dof*>(42);
}


#endif
