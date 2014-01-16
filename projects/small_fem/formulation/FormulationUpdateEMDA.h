#ifndef _FORMULATIONUPDATEEMDA_H_
#define _FORMULATIONUPDATEEMDA_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationUpdateEMDA
   @brief Update Formulation for FormulationEDMA

   Update Formulation for FormulationEDMA
 */

class FormulationUpdateEMDA: public Formulation<std::complex<double> >{
 private:
  // Wavenumber & Chi //
  double k;
  double chi;

  // Function Space & Basis //
  const FunctionSpaceScalar* fspace;
  const Basis*                basis;

  // Domain //
  const GroupOfElement* goe;

  // Quadrature //
  fullMatrix<double>* gC;
  fullVector<double>* gW;
  GroupOfJacobian*    jac;

  // DDM //
  const std::map<Dof, std::complex<double> >* solution;
  const std::map<Dof, std::complex<double> >* oldG;

 public:
  FormulationUpdateEMDA(const FunctionSpaceScalar& fs,
                        double k,
                        double chi,
                        const std::map<Dof, std::complex<double> >& solution,
                        const std::map<Dof, std::complex<double> >& oldG);

  virtual ~FormulationUpdateEMDA(void);

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
    interpolate(const MElement& element,
                const fullVector<double>& xyz,
                const std::map<Dof, std::complex<double> >& f) const;
};

/**
   @fn FormulationUpdateEMDA::FormulationUpdateEMDA
   @param goe A GroupOfElement
   @param k A real number
   @param chi A real number
   @param order A natural number
   @param ddmDof A map with the DDM Dof%s and their associated values

   Instantiates a new FormulationUpdateEMDA of the given order, wavenumber (k),
   real shift (chi) and ddm Dof%s

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationUpdateEMDA::~FormulationUpdateEMDA
   Deletes this FormulationUpdateEMDA
*/

#endif
