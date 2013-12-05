#ifndef _FORMULATIONPROJECTIONSCALAR_H_
#define _FORMULATIONPROJECTIONSCALAR_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "fullMatrix.h"

#include "TermFieldField.h"
#include "TermProjectionField.h"

#include "Formulation.h"

/**
   @class FormulationProjectionScalar
   @brief Formulation for the Projection of a Scalar Function problem

   Scalar Formulation for the L2 projection problem
 */

template<typename scalar>
class FormulationProjectionScalar: public Formulation<scalar>{
 private:
  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  const Basis*         basis;

  // For real version (Local Terms) //
  TermFieldField*      localTerms1;
  TermProjectionField* localTerms2;

  // For complex version //
  std::complex<double> (*f)(fullVector<double>& xyz);
  GroupOfElement*     goe;
  fullMatrix<double>* gC;
  fullVector<double>* gW;
  GroupOfJacobian*    jac;

 public:
  FormulationProjectionScalar(scalar (*f)(fullVector<double>& xyz),
                              FunctionSpaceScalar& fs);

  virtual ~FormulationProjectionScalar(void);

  virtual bool isGeneral(void) const;

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId)  const;
  virtual scalar weakB(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationProjectionScalar::FormulationProjectionScalar
   @param f The function to project
   @param fs A FunctionSpaceNode

   Instantiates a new FormulationProjectionScalar to project
   the given function@n

   FormulationProjectionScalar will use the given FunctionSpace
   for the projection
   **

   @fn FormulationProjectionScalar::~FormulationProjectionScalar
   Deletes the this FormulationProjectionScalar
*/

#endif
