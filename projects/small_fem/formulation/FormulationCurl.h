#ifndef _FORMULATIONCURL_H_
#define _FORMULATIONCURL_H_

#include "FunctionSpaceVector.h"
#include "TermCurlCurl.h"
#include "Formulation.h"

/**
   @class FormulationCurl
   @brief Formulation for the Curl problem

   Formulation for the @em Curl problem
 */

class FormulationCurl: public Formulation{
 private:
  // Function Space & Basis//
  FunctionSpaceVector* fspace;
  Basis*               basis;

  // Local Terms //
  TermCurlCurl* localTerms;

 public:
  FormulationCurl(GroupOfElement& goe, size_t order);

  virtual ~FormulationCurl(void);

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
   @fn FormulationCurl::FormulationCurl
   @param goe A GroupOfElement
   @param order A natural number

   Instantiates a new FormulationCurl of the given order@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationCurl::~FormulationCurl
   Deletes this FormulationCurl
*/

#endif
