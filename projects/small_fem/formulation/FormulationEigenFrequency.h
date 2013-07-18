#ifndef _FORMULATIONEIGENFREQUENCY_H_
#define _FORMULATIONEIGENFREQUENCY_H_

#include "FunctionSpaceVector.h"

#include "TermCurlCurl.h"
#include "TermGradGrad.h"

#include "Formulation.h"

/**
   @class FormulationEigenFrequency
   @brief Formulation for the Eigenfrequencies Problem

   Formulation for the Eigenfrequencies Problem
*/

class FormulationEigenFrequency: public Formulation{
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
                            size_t order);

  virtual ~FormulationEigenFrequency(void);

  virtual bool isGeneral(void) const;

  virtual double weak(size_t dofI, size_t dofJ,
                      size_t elementId) const;

  virtual double weakB(size_t dofI, size_t dofJ,
                       size_t elementId) const;

  virtual double rhs(size_t equationI,
                     size_t elementId) const;

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
