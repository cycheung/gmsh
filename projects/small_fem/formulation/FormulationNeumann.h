#ifndef _FORMULATIONNEUMANN_H_
#define _FORMULATIONNEUMANN_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationNeumann
   @brief Neumann Formulation

   Neumann Formulation
 */

class FormulationNeumann: public Formulation<std::complex<double> >{
 private:
  // Pulsation //
  double omega;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;

  // Local Terms //
  TermFieldField* localTerms;

 public:
  FormulationNeumann(GroupOfElement& goe,
                              double omega,
                              size_t order);

  virtual ~FormulationNeumann(void);

  virtual bool isGeneral(void) const;

  virtual std::complex<double>
    weak(size_t dofI, size_t dofJ, size_t elementId)  const;

  virtual std::complex<double>
    weakB(size_t dofI, size_t dofJ, size_t elementId) const;

  virtual std::complex<double>
    rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationNeumann::FormulationNeumann
   @param goe A GroupOfElement
   @param omega A real number
   @param order A natural number

   Instantiates a new FormulationNeumann of the given
   order and pulsation (omega)@n

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationNeumann::~FormulationNeumann
   Deletes this FormulationNeumann
*/

#endif
