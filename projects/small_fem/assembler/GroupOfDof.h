#ifndef _GROUPOFDOF_H_
#define _GROUPOFDOF_H_

#include <vector>
#include "Dof.h"
#include "Element.h"
#include "Jacobian.h"

/**
   @class GroupOfDof
   @brief Handel a group of Dof%s with @em geometrical meaning

   This class handles a group of Dof%s with a @em geometrical meaning 
   (@e e.g: Dof%s that belongs to the same (finite) element).@n

   Every GroupOfDof got a @em Jacobian transformaion matrix, and Dof%s @em orientations.@n
   These parameters are related to a geometrical element.@n
   The Dof%s of the GroupOfDof are also related to this element.@n
   For instance, the degrees of freedom of a finite element are a GroupOfDof.@n

   Moreover, every GroupOfDof has a unique (for a given DofManager) @c ID.

   @warning
   This class should be completly private to the @em DofManager and/or the @em System.@n
   The user shall not be aware of geometry, when it comes to Dof.

   @todo
   Chang GroupOfDof to a private class for the DofManager and/or System.
*/


class DofManager;

class GroupOfDof{
 private:
  int id;

  std::vector<Dof*>* dof;
  const std::vector<int>* direction;
  int nDof;
  int nextDof;
  
  Jacobian* jac;
  friend class DofManager;

 private:
   GroupOfDof(int numberOfDof, int groupId);
  ~GroupOfDof(void);

  void add(Dof* dof);
  void jacobian(Element& element);
  void orientation(const std::vector<int>& orientation);

 public:
  int getId(void) const;
  int dofNumber(void) const;
  const std::vector<Dof*>& getAllDofs(void) const;
  const Jacobian& getJacobian(void) const;
  int getOrientation(const int dofId) const;
};


/**
   @fn int GroupOfDof::getId(void) const
   @return Returns the @c ID of this GroupOfDof

   @fn int GroupOfDof::dofNumber(void) const
   @return Returns the number of Dof in this GroupOfDof

   @fn const std::vector<Dof*>& GroupOfDof::getAllDofs(void) const
   @return Returns all the Dof%s in this GroupOfDof

   @fn const Jacobian& GroupOfDof::getJacobian(void) const;
   @return Returns the Jacobian associated to this GroupOfDof
   
   @fn int GroupOfDof::getOrientation(const int dofId) const;
   @param dofId The @em local @c ID of a Dof in the GroupOfDof
   @return Returns the orientation of a Dof (indentified by its @em local @c ID)
*/

//////////////////////
// Inline Functions //
//////////////////////

inline void GroupOfDof::jacobian(Element& element){
  jac = &(element.getJacobian());
}

inline void GroupOfDof::orientation(const std::vector<int>& orientation){
  direction = &orientation;
}

inline int GroupOfDof::getId(void) const{
  return id;
}

inline int GroupOfDof::dofNumber(void) const{
  return nDof;
}

inline const std::vector<Dof*>& GroupOfDof::getAllDofs(void) const{
  return *dof;
}

inline const Jacobian& GroupOfDof::getJacobian(void) const{
  return *jac;
}

inline int GroupOfDof::getOrientation(const int dofId) const{
  return (*direction)[dofId];
}

#endif
