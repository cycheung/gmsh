#ifndef _FORMULATIONOO2_H_
#define _FORMULATIONOO2_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "Formulation.h"

/**
   @class FormulationOO2
   @brief OO2 Formulation for DDM

   OO2 Formulation for DDM
 */

class FormulationOO2: public Formulation<std::complex<double> >{
 private:
  // a & b //
  std::complex<double> a;
  std::complex<double> b;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*                basis;

  // Domain //
  GroupOfElement* goe;

  // Local Terms //
  TermFieldField* localTermsUU;
  TermGradGrad*   localTermsGG;

  // Quadrature (Field - Field) //
  fullMatrix<double>* gC;
  fullVector<double>* gW;
  GroupOfJacobian*    jac;

  // DDM //
  const std::map<Dof, std::complex<double> >* ddmDof;

 public:
  FormulationOO2(GroupOfElement& goe,
                 std::complex<double> a,
                 std::complex<double> b,
                 size_t order,
                 const std::map<Dof, std::complex<double> >& ddmDof);

  virtual ~FormulationOO2(void);

  virtual bool isGeneral(void) const;

  virtual std::complex<double>
    weak(size_t dofI, size_t dofJ, size_t elementId)  const;

  virtual std::complex<double>
    weakB(size_t dofI, size_t dofJ, size_t elementId) const;

  virtual std::complex<double>
    rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace& fs(void) const;

 private:
  std::complex<double>
    interpolate(const MElement& element, const fullVector<double>& xyz) const;
};

/**
   @fn FormulationOO2::FormulationOO2
   @param goe A GroupOfElement
   @param k A real number
   @param chi A real number
   @param order A natural number
   @param ddmDof A map with the DDM Dof%s and their associated values

   Instantiates a new FormulationOO2 of the given order,
   coefficients (a & b) and ddm Dof%s

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationOO2::~FormulationOO2
   Deletes this FormulationOO2
*/

#endif
