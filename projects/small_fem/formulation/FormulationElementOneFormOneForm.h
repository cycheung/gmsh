#ifndef _FORMULATIONELEMENTONEFORMONEFORM_H_
#define _FORMULATIONELEMENTONEFORMONEFORM_H_

#include "FunctionSpaceVector.h"
#include "Formulation.h"

/**
   @class FormulationElementOneFormOneForm
   @brief Formulation for an elementary 1-Form 1-Form matrix

   Formulation for an elementary 1-Form 1-Form matrix
 */

class FormulationElementOneFormOneForm: public Formulation{
 private:
  // Basis //
  FunctionSpaceVector* fspace;
  Basis*               basis;
  size_t               orientation;

  // Quadrature //
  fullMatrix<double>* gC;
  fullVector<double>* gW;

 public:
  FormulationElementOneFormOneForm(GroupOfElement& goe,
                                   size_t order,
                                   size_t orientation);

  virtual ~FormulationElementOneFormOneForm(void);

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
   @fn FormulationElementOneFormOneForm::FormulationElementOneFormOneForm
   @param geo A GroupOfElement defining the geomtry to use
   @param order The order of the basis to use
   @param orientation The orientation to use in the generated Basis

   Instantiates a new FormulationElementOneFormOneForm with the given options.
   **

   @fn FormulationElementOneFormOneForm::~FormulationElementOneFormOneForm
   Deletes this FormulationElementOneFormOneForm
*/

#endif
