#ifndef _FORMULATIONUPDATEOO2_H_
#define _FORMULATIONUPDATEOO2_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "Formulation.h"

/**
   @class FormulationUpdateOO2
   @brief Update Formulation for FormulationOO2

   Update Formulation for FormulationOO2
 */

class FormulationUpdateOO2: public Formulation<std::complex<double> >{
 private:
  // a & b //
  std::complex<double> a;
  std::complex<double> b;

  // Function Space & Basis //
  const FunctionSpaceScalar* fspace;
  const Basis*                basis;

  // Domain //
  const GroupOfElement* goe;

  // Quadrature (Field - Field) //
  fullMatrix<double>* gCFF;
  fullVector<double>* gWFF;
  GroupOfJacobian*    jacFF;

  // Quadrature (Grad - Grad) //
  fullMatrix<double>* gCGG;
  fullVector<double>* gWGG;
  GroupOfJacobian*    jacGG;

  // DDM //
  const std::map<Dof, std::complex<double> >* solution;
  const std::map<Dof, std::complex<double> >* oldG;

 public:
  FormulationUpdateOO2(const FunctionSpaceScalar& fs,
                       std::complex<double> a,
                       std::complex<double> b,
                       const std::map<Dof, std::complex<double> >& solution,
                       const std::map<Dof, std::complex<double> >& oldG);

  virtual ~FormulationUpdateOO2(void);

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

  fullVector<std::complex<double> >
    interpolateGrad(const MElement& element,
                    const fullVector<double>& xyz,
                    const std::map<Dof, std::complex<double> >& f) const;
};

/**
   @fn FormulationUpdateOO2::FormulationUpdateOO2
   @param goe A GroupOfElement
   @param a A real number
   @param b A real number
   @param order A natural number
   @param ddmDof A map with the DDM Dof%s and their associated values

   Instantiates a new FormulationUpdateOO2 of the given order,
   coefficients (a & b) and ddm Dof%s

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationUpdateOO2::~FormulationUpdateOO2
   Deletes this FormulationUpdateOO2
*/

#endif
