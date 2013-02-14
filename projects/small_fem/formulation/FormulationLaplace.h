#ifndef _FORMULATIONLAPLACE_H_
#define _FORMULATIONLAPLACE_H_

#include "FunctionSpaceScalar.h"

#include "TermGradGrad.h"

#include "Formulation.h"

/**
   @class FormulationLaplace
   @brief Formulation for the Laplace problem

   Formulation for the @em Laplace problem
 */

class FormulationLaplace: public Formulation{
 private:
  // Function Space & Basis//
  FunctionSpaceScalar* fspace;
  Basis*               basis;

  // Local Terms //
  TermGradGrad* localTerms;

 public:
  FormulationLaplace(GroupOfElement& goe, unsigned int order);

  virtual ~FormulationLaplace(void);

  virtual bool isGeneral(void) const;

  virtual double weak(unsigned int dofI, unsigned int dofJ,
		      const GroupOfDof& god) const;

  virtual double weakB(unsigned int dofI, unsigned int dofJ,
                       const GroupOfDof& god) const;

  virtual double rhs(unsigned int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationLaplace::FormulationLaplace
   @param goe A GroupOfElement
   @param order A natural number

   Instantiates a new FormulationLaplace of the given order@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationLaplace::~FormulationLaplace
   Deletes this FormulationLaplace
   **
*/

#endif
