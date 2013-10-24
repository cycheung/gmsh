#ifndef _FORMULATIONELEMENTGRADGRAD_H_
#define _FORMULATIONELEMENTGRADGRAD_H_

#include "FunctionSpaceScalar.h"
#include "Formulation.h"

/**
   @class FormulationElementGradGrad
   @brief Formulation for an elementary grad grad matrix

   Formulation for an elementary grad grad matrix
 */

class FormulationElementGradGrad: public Formulation{
 private:
  // Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;
  size_t               orientation;

  // Quadrature //
  fullMatrix<double>* gC;
  fullVector<double>* gW;

 public:
  FormulationElementGradGrad(GroupOfElement& goe,
                             size_t order,
                             size_t orientation);

  virtual ~FormulationElementGradGrad(void);

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
   @fn FormulationElementGradGrad::FormulationElementGradGrad
   @param geo A GroupOfElement defining the geomtry to use
   @param order The order of the basis to use
   @param orientation The orientation to use in the generated Basis

   Instantiates a new FormulationElementGradGrad with the given options.
   **

   @fn FormulationElementGradGrad::~FormulationElementGradGrad
   Deletes this FormulationElementGradGrad
*/

#endif
