#ifndef _FORMULATIONPROJECTIONVECTOR_H_
#define _FORMULATIONPROJECTIONVECTOR_H_

#include "FunctionSpaceVector.h"
#include "fullMatrix.h"

#include "TermGradGrad.h"
#include "TermProjectionGrad.h"

#include "Formulation.h"

/**
   @class FormulationProjectionVector
   @brief Formulation for the Projection of Vectorial Function problem

   Vectorial Formulation for the @em L2 @em Projection problem
 */

class FormulationProjectionVector: public Formulation{
 private:
  // Function to Project //
  fullVector<double> (*f)(fullVector<double>& xyz);

  // Function Space & Basis //
  FunctionSpaceVector* fspace;
  const Basis*         basis;

  // Local Terms //
  TermGradGrad*       localTerms1;
  TermProjectionGrad* localTerms2;

 public:
  FormulationProjectionVector(fullVector<double> (*f)(fullVector<double>& xyz),
			      FunctionSpaceVector& fs);

  virtual ~FormulationProjectionVector(void);

  virtual bool isGeneral(void) const;

  virtual double weak(unsigned int dofI, unsigned int dofJ,
		      unsigned int elementId) const;

  virtual double weakB(unsigned int dofI, unsigned int dofJ,
                       unsigned int elementId) const;

  virtual double rhs(unsigned int equationI,
		     unsigned int elementId) const;

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
