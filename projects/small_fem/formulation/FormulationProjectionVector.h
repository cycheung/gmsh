#ifndef _FORMULATIONPROJECTIONVECTOR_H_
#define _FORMULATIONPROJECTIONVECTOR_H_

#include "FunctionSpaceVector.h"
#include "fullMatrix.h"

#include "TermHCurl.h"
#include "TermProjectionHCurl.h"

#include "Formulation.h"

/**
   @class FormulationProjectionVector
   @brief Formulation for the Projection of Vectorial Function problem

   Vectorial Formulation for the @em L2 @em Projection problem.
 */

class FormulationProjectionVector: public Formulation{
 private:
  // Function to Project //
  fullVector<double> (*f)(fullVector<double>& xyz);

  // Function Space & Basis //
  FunctionSpaceVector* fspace;
  const Basis*         basis;

  // Local Terms //
  TermHCurl*           localTerms1;
  TermProjectionHCurl* localTerms2;

 public:
  FormulationProjectionVector(fullVector<double> (*f)(fullVector<double>& xyz),
			      FunctionSpaceVector& fs);

  virtual ~FormulationProjectionVector(void);

  virtual double weak(unsigned int dofI, unsigned int dofJ,
		      const GroupOfDof& god) const;

  virtual double rhs(unsigned int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationProjectionVector::FormulationProjectionVector
   @param f The function to project
   @param fs A FunctionSpaceEdge

   Instantiates a new FormulationProjectionVector to project
   the given function@n

   FormulationProjectionVector will use the given FunctionSpace
   for the projection
   **

   @fn FormulationProjectionVector::~FormulationProjectionVector
   Deletes the this FormulationProjectionVector
   **
*/

#endif
