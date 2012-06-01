#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include "fullMatrix.h"
#include "Mesh.h"
#include "Jacobian.h"
#include "DofManager.h"
#include "Formulation.h"

/**
   @class System
   @brief This class assembles a linear system

   This class assembles the linear system, that comes from a Formulation 
   and a list of Element%s.@n
  
   @warning
   Up to now, the assembly is done when the System is instantiated@n
   Also, we can @em only assemble Dof related to an Element@n

   @todo
   Assembly done by a specific method, not the constructor@n
   Assembly of @em non @em geometric Dof@n
   Maybe put the list of element on the formulation ? 
   Or more abstract concept -- 'Problem' Class ?@n
   Give the possiblite to get the solution Vector instead of setting 
   Entity values
 */

class System{
 private:
  fullMatrix<double>* A;
  fullVector<double>* n;
  int size;

  DofManager* dofM;
  const Formulation* formulation;

 public:
   System(const std::vector<Element*>& elements, 
	  const Formulation& formulation);
  ~System(void);

  fullMatrix<double>& getMatrix(void) const;
  fullVector<double>& getRHS(void) const;

  void fixBC(const int physicalId, const double value);
  void solve(void);
  
 private:
  void assemble(GroupOfDof& group);
};


/**
   @fn System::System(const std::vector<Element*>& elements, 
   const Formulation& formulation)
   @param elements A list of Element%s, giving the geomtry of the problem to solve
   @param formulation A Formulation, giving the way to assemble the system
   @return A new @em assembled System
 
   @fn System::~System(void)
   @return Deletes the System

   @fn Matrix& System::getMatrix(void) const
   @return Returns the assembled matrix of the the linear system

   @fn fullVector<double>& System::getRHS(void) const
   @return Returns the assembled Right Hand Side of the the linear system

   @fn void System::fixBC(const int physicalId, const double value)
   @param physicalId The physical @c ID on which the bondary condtion shall be imposed
   @param value The value of the bondary condition
   @return Fix a Boundary Condition on the linear system

   @fn void System::solve(void)
   @return Solves the linear system
*/

//////////////////////
// Inline Functions //
//////////////////////


inline fullMatrix<double>& System::getMatrix(void) const{
  return *A;
}

inline fullVector<double>& System::getRHS(void) const{
  return *n;
}

#endif
