#ifndef _FORMULATIONSTEADYWAVEVECTORSLOW_H_
#define _FORMULATIONSTEADYWAVEVECTORSLOW_H_

#include <vector>

#include "FunctionSpaceVector.h"
#include "GroupOfJacobian.h"
#include "Formulation.h"

/**
   @class FormulationSteadyWaveVectorSlow
   @brief Vectorial Formulation for the Steady Wave problem (Slow version)

   Same as FormulationSteadyWaveVector, but this version don't use
   the fast integration algorithm
 */

class FormulationSteadyWaveVectorSlow: public Formulation{
 private:
  // Physical Values //
  static const double mu;
  static const double eps;

  // Wave Number Squared //
  double kSquare;

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
                                  double k,
                                  unsigned int order);

  virtual ~FormulationSteadyWaveVectorSlow(void);


  virtual bool isGeneral(void) const;

  virtual double weak(unsigned int dofI, unsigned int dofJ,
		      unsigned int elementId) const;

  virtual double weakB(unsigned int dofI, unsigned int dofJ,
                       unsigned int elementId) const;

  virtual double rhs(unsigned int equationI,
		     unsigned int elementId) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationSteadyWaveVectorSlow::FormulationSteadyWaveVectorSlow
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWaveVectorSlow of the given
   @em order and @em wave @em number (@c k)@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationSteadyWaveVectorSlow::~FormulationSteadyWaveVectorSlow
   Deletes this FormulationSteadyWaveVectorSlow
*/

#endif
