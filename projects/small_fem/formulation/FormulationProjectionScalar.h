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
  FormulationProjectionScalar(double (*f)(fullVector<double>& xyz),
			      FunctionSpaceNode& fs);
  
  virtual ~FormulationProjectionScalar(void);

  virtual double weak(int dofI, int dofJ, 
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationProjectionScalar::FormulationProjectionScalar
   @param f The function to project
   @param fs A FunctionSpaceNode
   
   Instantiates a new FormulationProjectionScalar to project
   the given function@n

   FormulationProjectionScalar will use the given FunctionSpace
   for the projection

   @warning The given FunctionSpace will be ask to 
   @em PreEvaluate 
   (@em see FunctionSpaceScalar::preEvaluateLocalFunctions() 
   and similar)
   some Functions
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
