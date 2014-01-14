#ifndef _FORMULATIONEMDA_H_
#define _FORMULATIONEMDA_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "Formulation.h"

/**
   @class FormulationEMDA
   @brief EMDA Formulation for DDM

   EMDA Formulation for DDM
 */

class FormulationEMDA: public Formulation<std::complex<double> >{
 private:
  // Function Space & Basis //
  const FunctionSpaceScalar* fspace;
  const Basis*                basis;

  // Domain //
  const GroupOfElement* goe;

  // Qudrature //
  fullMatrix<double>* gC;
  fullVector<double>* gW;
  GroupOfJacobian*    jac;

  // DDM //
  const std::map<Dof, std::complex<double> >* ddmDof;

 public:
  FormulationEMDA(const FunctionSpaceScalar& fs,
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
   @param fs A FunctionSpaceScalar
   @param ddmDof A map with the DDM Dof%s and their associated values

   Instantiates a new FormulationEMDA with the given FunctionSpace and ddm Dof%s
   **

   @fn FormulationEMDA::~FormulationEMDA
   Deletes this FormulationEMDA
*/

#endif
