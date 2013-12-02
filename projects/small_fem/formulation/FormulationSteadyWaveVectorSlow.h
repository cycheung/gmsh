#ifndef _FORMULATIONSTEADYWAVEVECTORSLOW_H_
#define _FORMULATIONSTEADYWAVEVECTORSLOW_H_

#include <vector>

#include "FunctionSpaceVector.h"
#include "GroupOfJacobian.h"
#include "FormulationTyped.h"

/**
   @class FormulationSteadyWaveVectorSlow
   @brief Vectorial Formulation for the Steady Wave problem (Slow version)

   Same as FormulationSteadyWaveVector, but this version don't use
   the fast integration algorithm
 */

class FormulationSteadyWaveVectorSlow: public FormulationTyped<double>{
 private:
  // Speed of medium squared //
  static const double cSquare;

  // Pulsation Squared //
  double omegaSquare;

  // Gaussian Quadrature Data (Term One) //
  int G1;
  fullMatrix<double>* gC1;
  fullVector<double>* gW1;

  // Gaussian Quadrature Data (Term Two) //
  int G2;
  fullMatrix<double>* gC2;
  fullVector<double>* gW2;

  // Domain //
  GroupOfElement* goe;

  // Jacobians //
  GroupOfJacobian* jac1;
  GroupOfJacobian* jac2;

  // Function Space & Basis //
  FunctionSpaceVector* fspace;
  Basis*               basis;

 public:
  FormulationSteadyWaveVectorSlow(GroupOfElement& goe,
                                  double omega,
                                  size_t order);

  virtual ~FormulationSteadyWaveVectorSlow(void);


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
   @fn FormulationSteadyWaveVectorSlow::FormulationSteadyWaveVectorSlow
   @param goe A GroupOfElement
   @param omega A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWaveVectorSlow of the given
   order and pulsation (omega)@n

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationSteadyWaveVectorSlow::~FormulationSteadyWaveVectorSlow
   Deletes this FormulationSteadyWaveVectorSlow
*/

#endif
