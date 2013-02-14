#ifndef _SLOWFORMULATIONSTEADYWAVEVECTOR_H_
#define _SLOWFORMULATIONSTEADYWAVEVECTOR_H_

#include "FunctionSpaceVector.h"
#include "Jacobian.h"
#include "Formulation.h"

/**
   @class SlowFormulationSteadyWaveVector
   @brief A Slower Version of FormulationSteadyWaveVector

   A Slower Version of FormulationSteadyWaveVector

   @see FormulationSteadyWaveVector
 */

class SlowFormulationSteadyWaveVector: public Formulation{
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

  // Function Space & Basis //
  FunctionSpaceVector* fspace;
  Basis*               basis;

  // Jacobians //
  Jacobian* jac1;
  Jacobian* jac2;

 public:
  SlowFormulationSteadyWaveVector(GroupOfElement& goe,
                                  double k,
                                  unsigned int order);

  virtual ~SlowFormulationSteadyWaveVector(void);

  virtual bool isGeneral(void) const;

  virtual double weak(unsigned int dofI, unsigned int dofJ,
		      const GroupOfDof& god) const;

  virtual double weakB(unsigned int dofI, unsigned int dofJ,
                       const GroupOfDof& god) const;

  virtual double rhs(unsigned int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};


#endif
