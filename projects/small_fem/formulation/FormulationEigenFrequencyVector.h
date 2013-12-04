#ifndef _FORMULATIONEIGENFREQUENCYVECTOR_H_
#define _FORMULATIONEIGENFREQUENCYVECTOR_H_

#include <complex>
#include "FunctionSpaceVector.h"

#include "TermCurlCurl.h"
#include "TermGradGrad.h"

#include "Formulation.h"

/**
   @class FormulationEigenFrequencyVector
   @brief Formulation for the vectrorial Eigenfrequencies Problem

   Formulation for the vectorial Eigenfrequencies Problem
*/

class FormulationEigenFrequencyVector:
public Formulation<std::complex<double> >{
 private:
  // Speed of medium squared //
  static const double cSquare;

  // Function Space & Basis //
  FunctionSpaceVector* fspace;
  Basis*               basis;

  // Local Terms //
  TermCurlCurl* localTerms1;
  TermGradGrad* localTerms2;

 public:
  FormulationEigenFrequencyVector(GroupOfElement& goe,
                            size_t order);

  virtual ~FormulationEigenFrequencyVector(void);

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
   @fn FormulationEigenFrequencyVector::FormulationEigenFrequencyVector
   @param goe A GroupOfElement defining the Domain of the Problem
   @param order A natural number, giving the order of this Formulation

   Instanciates a new Formulation for the vectorial
   Eigenfrequencies Problem
   **

   @fn FormulationEigenFrequencyVector::~FormulationEigenFrequencyVector
   Deletes this Formualtion
*/

#endif
