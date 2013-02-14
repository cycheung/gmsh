#ifndef _FORMULATIONPOISSON_H_
#define _FORMULATIONPOISSON_H_

#include "FunctionSpaceScalar.h"
#include "fullMatrix.h"

#include "TermGradGrad.h"
#include "TermProjectionField.h"

#include "Formulation.h"

/**
   @class FormulationPoisson
   @brief Formulation for the Poisson problem

   Formulation for the @em Poisson problem
 */

class FormulationPoisson: public Formulation{
 private:
  // Source Term //
  static const unsigned int sourceOrder;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;

  // Local Terms //
  TermGradGrad*        localTermsL;
  TermProjectionField* localTermsR;

 public:
  FormulationPoisson(GroupOfElement& goe,
		     unsigned int order);

  virtual ~FormulationPoisson(void);

  virtual bool isGeneral(void) const;

  virtual double weak(unsigned int dofI, unsigned int dofJ,
		      unsigned int elementId) const;

  virtual double weakB(unsigned int dofI, unsigned int dofJ,
                       unsigned int elementId) const;

  virtual double rhs(unsigned int equationI,
		     unsigned int elementId) const;

  virtual const FunctionSpace& fs(void) const;

 private:
  static double gSource(fullVector<double>& xyz);
};

/**
   @fn FormulationPoisson::FormulationPoisson
   @param goe A GroupOfElement
   @param order A natural number

   Instantiates a new FormulationPoisson of the given order@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationPoisson::~FormulationPoisson
   Deletes this FormulationPoisson
*/

#endif
