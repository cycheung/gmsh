#ifndef _GROUPOFDOF_H_
#define _GROUPOFDOF_H_

#include <vector>
#include "Dof.h"
#include "Mapper.h"
#include "MElement.h"
#include "fullMatrix.h"

/**
   @class GroupOfDof
   @brief Handel a group of Dof%s with @em geometrical meaning

   This class handles a group of Dof%s with a @em geometrical meaning 
   (@e e.g: Dof%s that belongs to the same (finite) element).@n

   This class gives acces to individual Dof%s of the group.@n
   It also gives acces to the underlying Geometrical Element.  
*/


class DofManager;

class GroupOfDof{
 private:
  const MElement* element;

  int nDof;
  std::vector<Dof*>* dof;
  const std::vector<int>* direction;
  
  int nextDof;
  
  friend class DofManager;

 public:
  int getId(void) const;
  int dofNumber(void) const;
  const std::vector<Dof*>& getAllDofs(void) const;
  const MElement&          getElement(void) const;

  int getOrientation(const int dofId) const;

 private:
   GroupOfDof(int numberOfDof, const MElement& geoElement);
  ~GroupOfDof(void);

  void add(Dof* dof);
  void orientation(const std::vector<int>& orientation);
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

inline void GroupOfDof::orientation(const std::vector<int>& orientation){
  direction = &orientation;
}

inline int GroupOfDof::getId(void) const{
  return element->getNum();
}

inline int GroupOfDof::dofNumber(void) const{
  return nDof;
}

inline const std::vector<Dof*>& GroupOfDof::getAllDofs(void) const{
  return *dof;
}

inline const MElement& GroupOfDof::getElement(void) const{
  return *element;
}

inline int GroupOfDof::getOrientation(const int dofId) const{
  return (*direction)[dofId];
}

#endif
