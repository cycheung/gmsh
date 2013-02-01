#ifndef _FORMULATIONEIGENFREQUENCY_H_
#define _FORMULATIONEIGENFREQUENCY_H_

#include "FunctionSpaceVector.h"

#include "TermCurlCurl.h"
#include "TermGradGrad.h"

#include "EigenFormulation.h"

/**
   @class FormulationEigenFrequency
   @brief EigenFormulation for the Eigenfrequencies Problem

   EigenFormulation for the Eigenfrequencies Problem
 */

class FormulationEigenFrequency: public EigenFormulation{
 private:
  // Physical Values //
  static const double mu;
  static const double eps;

  // Function Space & Basis //
  FunctionSpaceVector* fspace;
  Basis*               basis;

  // Local Terms //
  TermCurlCurl* localTerms1;
  TermGradGrad* localTerms2;

 public:
  FormulationEigenFrequency(GroupOfElement& goe,
			    unsigned int order);

  virtual ~FormulationEigenFrequency(void);

  virtual bool isGeneral(void) const;

  virtual double weakA(unsigned int dofI, unsigned int dofJ,
		       const GroupOfDof& god) const;

  virtual double weakB(unsigned int dofI, unsigned int dofJ,
		       const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationEigenFrequency::FormulationEigenFrequency
   @param goe A GroupOfElement defining the Domain of the Problem
   @param order A natural number, giving the order of this Formulation

   Instanciates a new EigenFormulation for the
   Eigenfrequencies Problem
   **

   @fn FormulationEigenFrequency::~FormulationEigenFrequency
   Deletes this FormualtionEigenFrequency
*/

#endif
