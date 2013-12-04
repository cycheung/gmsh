#ifndef _FORMULATIONEIGENFREQUENCYSCALAR_H_
#define _FORMULATIONEIGENFREQUENCYSCALAR_H_

#include <complex>
#include "FunctionSpaceScalar.h"

#include "TermGradGrad.h"
#include "TermFieldField.h"

#include "Formulation.h"

/**
   @class FormulationEigenFrequencyScalar
   @brief Formulation for the scalar Eigenfrequencies Problem

   Formulation for the scalar Eigenfrequencies Problem
 */

class FormulationEigenFrequencyScalar:
public Formulation<std::complex<double> >{
 private:
  // Speed of medium squared //
  static const double cSquare;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;

  // Local Terms //
  TermGradGrad*   localTerms1;
  TermFieldField* localTerms2;

 public:
  FormulationEigenFrequencyScalar(GroupOfElement& goe,
                       size_t order);

  virtual ~FormulationEigenFrequencyScalar(void);

  virtual bool isGeneral(void) const;

  virtual std::complex<double> weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const;

  virtual std::complex<double> weakB(size_t dofI, size_t dofJ,
                                     size_t elementId) const;

  virtual std::complex<double> rhs(size_t equationI,
                                   size_t elementId) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationEigenFrequencyScalar::FormulationEigenFrequencyScalar
   @param goe A GroupOfElement
   @param order A natural number

   Instantiates a new FormulationEigenFrequencyScalar of the given order@n

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationEigenFrequencyScalar::~FormulationEigenFrequencyScalar
   Deletes this FormulationEigenFrequencyScalar
*/

#endif
