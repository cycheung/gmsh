#ifndef _FORMULATIONEMDA_H_
#define _FORMULATIONEMDA_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationEMDA
   @brief EMDA Formulation for DDM

   EMDA Formulation for DDM
 */

class FormulationEMDA: public Formulation<std::complex<double> >{
 private:
  // Wavenumber //
  double k;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*                basis;

  // Domain //
  GroupOfElement* goe;

  // Local Terms (Projection u_i*u_j) //
  TermFieldField* localTerms;

  // Quadrature //
  fullMatrix<double>* gC;
  fullVector<double>* gW;
  GroupOfJacobian*    jac;

  // DDM //
  const std::map<Dof, std::complex<double> >* ddmDof;

 public:
  FormulationEMDA(GroupOfElement& goe,
                  double k,
                  size_t order,
                  const std::map<Dof, std::complex<double> >& ddmDof);

  virtual ~FormulationEMDA(void);

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
   @fn FormulationEMDA::FormulationEMDA
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number
   @param ddmDof A map with the DDM Dof%s and their associated values

   Instantiates a new FormulationEMDA of the given order, wavenumber (k)
   and ddm Dof%s

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationEMDA::~FormulationEMDA
   Deletes this FormulationEMDA
*/

#endif
