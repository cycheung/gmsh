#ifndef _FORMULATIONPROJECTION_H_
#define _FORMULATIONPROJECTION_H_

#include <vector>

#include "FunctionSpaceEdge.h"
#include "fullMatrix.h"
#include "Formulation.h"

/**
   @class FormulationProjection
   @brief Formulation for the Projection problem

   Vectorial Formulation for the @em L2 @em Projection problem.
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
  FunctionSpaceEdge* fspace;
  
 public:
  FormulationProjection(const GroupOfElement& goe,
			const fullVector<double>& vectorToProject);
  
  virtual ~FormulationProjection(void);

  virtual double weak(int edgeI, int edgeJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationProjection::FormulationProjection
   @param goe A GroupOfElement
   @param vectorToProject A fullVector<double>
   
   Instantiates a new FormulationProjection to project
   the given Vector@n

   The given GroupOfElement will be used as the 
   geomtrical @em domain
   **

   @fn FormulationProjection::~FormulationProjection
   Deletes the this FormulationProjection
   **
*/


//////////////////////
// Inline Functions //
//////////////////////

inline const FunctionSpace& FormulationProjection::fs(void) const{
  return *fspace;
}

#endif
