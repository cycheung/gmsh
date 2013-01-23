#ifndef _FORMULATIONPOISSON_H_
#define _FORMULATIONPOISSON_H_

#include <vector>

#include "FunctionSpaceScalar.h"
#include "Formulation.h"

/**
   @class FormulationPoisson
   @brief Formulation for the Poisson problem

   Formulation for the @em Poisson problem.

   @todo
   Remove ALL const_cast%S by correcting MElement constness@n
   Allow Hybrid Mesh
 */

class FormulationPoisson: public Formulation{
 private:
  // Gaussian Quadrature Data (LHS) //
  int GL;
  fullMatrix<double>* gCL;
  fullVector<double>* gWL;

  // Gaussian Quadrature Data (RHS) //
  int GR;
  fullMatrix<double>* gCR;
  fullVector<double>* gWR;

  // Function Space & Basis //
  FunctionSpaceScalar* fspace;
  Basis*               basis;

 public:
  FormulationPoisson(GroupOfElement& goe,
		     unsigned int order);

  virtual ~FormulationPoisson(void);

  virtual double weak(int nodeI, int nodeJ,
		      const GroupOfDof& god) const;

  virtual double rhs(int equationI,
		     const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationPoisson::FormulationPoisson
   @param goe A GroupOfElement
   @param order A natural number

   Instantiates a new FormulationPoisson of the given order@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationPoisson::~FormulationPoisson
   Deletes this FormulationPoisson
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline const FunctionSpace& FormulationPoisson::fs(void) const{
  return *fspace;
}

#endif
