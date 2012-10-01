#ifndef _FORMULATIONPROJECTIONSCALAR_H_
#define _FORMULATIONPROJECTIONSCALAR_H_

#include "FunctionSpaceNode.h"
#include "fullMatrix.h"
#include "Formulation.h"

/**
   @class FormulationProjectionScalar
   @brief Formulation for the Projection of a Scalar Function problem

   Scalar Formulation for the @em L2 @em Projection problem.
 */

class FormulationProjectionScalar: public Formulation{
 private:
  // Gaussian Quadrature Data //
  int G;
  fullMatrix<double>* gC;
  fullVector<double>* gW;

  // Function Space //
  FunctionSpaceNode* fspace;

  // Function to Project //
  double (*f)(fullVector<double>& xyz);

 public:
  FormulationProjectionScalar(const GroupOfElement& goe,
			      double (*f)(fullVector<double>& xyz),
			      unsigned int ordre);
  
  virtual ~FormulationProjectionScalar(void);

  virtual double weak(int dofI, int dofJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationProjectionScalar::FormulationProjectionScalar
   @param goe A GroupOfElement
   @param f The function to project
   @param ordre A strictly positive natural number
   
   Instantiates a new FormulationProjectionScalar to project
   the given function@n

   FormulationProjectionScalar will use a @em Nodal Basis
   of the given ordre@n

   The given GroupOfElement will be used as the 
   geomtrical @em domain
   **

   @fn FormulationProjectionScalar::~FormulationProjectionScalar
   Deletes the this FormulationProjectionScalar
   **
*/


//////////////////////
// Inline Functions //
//////////////////////

inline const FunctionSpace& FormulationProjectionScalar::fs(void) const{
  return *fspace;
}

#endif
