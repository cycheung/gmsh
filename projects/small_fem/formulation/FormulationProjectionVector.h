#ifndef _FORMULATIONPROJECTIONVECTOR_H_
#define _FORMULATIONPROJECTIONVECTOR_H_

#include <vector>

#include "FunctionSpaceEdge.h"
#include "fullMatrix.h"
#include "Formulation.h"

/**
   @class FormulationProjectionVector
   @brief Formulation for the Projection of Vectorial Function problem

   Vectorial Formulation for the @em L2 @em Projection problem.
 */

class FormulationProjectionVector: public Formulation{
 private:
  // Gaussian Quadrature Data //
  int G;
  fullMatrix<double>* gC;
  fullVector<double>* gW;

  // Function Space //
  FunctionSpaceEdge* fspace;

  // Function to Project //
  fullVector<double> (*f)(fullVector<double>& xyz);
  
 public:
  FormulationProjectionVector(const GroupOfElement& goe,
			      fullVector<double> (*f)(fullVector<double>& xyz),
			      unsigned int order);
  
  virtual ~FormulationProjectionVector(void);

  virtual double weak(int dofI, int dofJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationProjectionVector::FormulationProjectionVector
   @param goe A GroupOfElement
   @param vectorToProject A fullVector<double>
   
   Instantiates a new FormulationProjectionVector to project
   the given Vector@n

   The given GroupOfElement will be used as the 
   geomtrical @em domain
   **

   @fn FormulationProjectionVector::~FormulationProjectionVector
   Deletes the this FormulationProjectionVector
   **
*/


//////////////////////
// Inline Functions //
//////////////////////

inline const FunctionSpace& FormulationProjectionVector::fs(void) const{
  return *fspace;
}

#endif
