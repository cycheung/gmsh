#ifndef _FORMULATIONELEMENTCURLCURL_H_
#define _FORMULATIONELEMENTCURLCURL_H_

#include "FunctionSpaceVector.h"
#include "Formulation.h"

/**
   @class FormulationElementCurlCurl
   @brief Formulation for an elementary grad grad matrix

   Formulation for an elementary curl curl matrix
 */

class FormulationElementCurlCurl: public Formulation{
 private:
  // Basis //
  FunctionSpaceVector* fspace;
  Basis*               basis;
  size_t               orientation;

  // Quadrature //
  fullMatrix<double>* gC;
  fullVector<double>* gW;

 public:
    FormulationElementCurlCurl(GroupOfElement& goe,
                               size_t order,
                               size_t orientation);

  virtual ~FormulationElementCurlCurl(void);

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
   @fn FormulationElementCurlCurl::FormulationElementCurlCurl
   @param geo A GroupOfElement defining the geomtry to use
   @param order The order of the basis to use
   @param orientation The orientation to use in the generated Basis

   Instantiates a new FormulationElementCurlCurl with the given options.
   **

   @fn FormulationElementCurlCurl::~FormulationElementCurlCurl
   Deletes this FormulationElementCurlCurl
*/

#endif
