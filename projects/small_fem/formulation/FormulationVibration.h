#ifndef _FORMULATIONVIBRATION_H_
#define _FORMULATIONVIBRATION_H_

#include <vector>

#include "FunctionSpaceNode.h"
#include "EigenFormulation.h"

/**
   @class FormulationVibration
   @brief EigenFormulation for the Vibration problem

   EigenFormulation for the @em Vibration problem.

   @todo
   Remove ALL const_cast%S by correcting MElement constness@n
   Allow Hybrid Mesh
 */

class FormulationVibration: public EigenFormulation{
 private:
  // Gaussian Quadrature Data (LHS) //
  int GL;
  fullMatrix<double>* gCL;
  fullVector<double>* gWL;

  // Function Space & Basis //
  FunctionSpaceNode* fspace;
  const Basis*       basis;

 public:
  FormulationVibration(const GroupOfElement& goe,
		       unsigned int order);

  virtual ~FormulationVibration(void);

  virtual bool isGeneral(void) const;

  virtual double weakA(int dofI, int dofJ,
		       const GroupOfDof& god) const;

  virtual double weakB(int dofI, int dofJ,
		       const GroupOfDof& god) const;

  virtual const FunctionSpace& fs(void) const;
};

/**
   @fn FormulationVibration::FormulationVibration
   @param goe A GroupOfElement
   @param order A natural number

   Instantiates a new FormulationVibration of the given order@n

   The given GroupOfElement will be used as the
   geomtrical @em domain
   **

   @fn FormulationVibration::~FormulationVibration
   Deletes this FormulationVibration
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline bool FormulationVibration::isGeneral(void) const{
  return false;
}

inline double FormulationVibration::weakB(int dofI, int dofJ,
					  const GroupOfDof& god) const{
  return 42;
}

inline const FunctionSpace& FormulationVibration::fs(void) const{
  return *fspace;
}

#endif
