#ifndef _FORMULATIONPROJECTION_H_
#define _FORMULATIONPROJECTION_H_

#include <vector>

#include "Formulation.h"
#include "Polynomial.h"
#include "fullMatrix.h"

/**
   @class FormulationProjection
   @brief Formulation for the Projection problem

   Vectorial Formulation for the @em Projection problem.
 */

class FormulationProjection: public Formulation{
 private:
  // Vector to Project //
  const fullVector<double>* f;

  // Gaussian Quadrature Data //
  int G;
  fullMatrix<double>* gC;
  fullVector<double>* gW;

  // Function Space //
  FunctionSpace*                               fspace;
  const std::vector<std::vector<Polynomial> >* basis;

 public:
  FormulationProjection(const GroupOfElement& goe,
			const fullVector<double>& vectorToProject);
  
  virtual ~FormulationProjection(void);

  virtual double weak(int edgeI, int edgeJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual FunctionSpace& fs(void) const;
};

/**
   @fn FormulationProjection::FormulationProjection
   @param vectorToProject A fullVector<double>
   @return Returns a new FormulationProjection to project
   the given Vector
 
   @fn FormulationProjection::~FormulationProjection
   @return Deletes the this FormulationProjection
*/


//////////////////////
// Inline Functions //
//////////////////////

inline FunctionSpace& FormulationProjection::fs(void) const{
  return *fspace;
}

#endif
