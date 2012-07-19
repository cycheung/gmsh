#ifndef _GROUPOFDOF_H_
#define _GROUPOFDOF_H_

#include <vector>
#include "Group.h"
#include "Dof.h"
#include "Mapper.h"
#include "MElement.h"
#include "fullMatrix.h"

class DofManager;

/**
   @class GroupOfDof
   @brief Handels a Group of Dof%s with @em geometrical meaning

   This class handles a group of Dof%s with a @em geometrical meaning 
   (@e e.g: Dof%s that belongs to the same (finite) element).@n

   This class gives acces to individual Dof%s of the group.@n
   It also gives acces to the underlying Geometrical Element.@n

   To conclude, this class is a Group.

   @note
   Note that a user @em can't instantiate a GroupOfDof.@n
   This is the DofManager job.
*/

class GroupOfDof: public Group{
 private:
  const MElement* element;

  int nDof;
  std::vector<Dof*>* dof;
  const std::vector<int>* direction;
  
  int nextDof;
  
  friend class DofManager;

 public:
  virtual int getNumber(void) const;
  virtual int getId(void)     const;
  virtual int getType(void)   const;
  
  Dof&                     get(int i) const; 
  const std::vector<Dof*>& getAll(void) const;
  
  const MElement&          getGeoElement(void) const;

  int getOrientation(const int dofId) const;

  virtual std::string toString(void) const;

 private:
   GroupOfDof(int numberOfDof, const MElement& geoElement);
  ~GroupOfDof(void);

  void add(Dof* dof);
  void orientation(const std::vector<int>& orientation);
};


/** 
   @fn GroupOfDof::get
   @param i An interger ranging from 0 
   to GroupOfDof::getNumber() - 1
   @return Returns the ith element of the Group

   @fn GroupOfDof::getAll
   @return Returns all the elements of the Group

   @fn int GroupOfDof::getOrientation(const int dofId) const;
   @param dofId The @em local @c ID of a Dof in the GroupOfDof
   @return Returns the orientation of a Dof (indentified by its @em local @c ID)
*/

//////////////////////
// Inline Functions //
//////////////////////

inline int GroupOfDof::getNumber(void) const{
  return nDof;
}

inline int GroupOfDof::getId(void) const{
  return element->getNum();
}

inline int GroupOfDof::getType(void) const{
  return 2;
}

inline Dof& GroupOfDof::get(int i) const{
  return *((*dof)[i]);
}
 
inline const std::vector<Dof*>& GroupOfDof::getAll(void) const{
  return *dof;
}

inline const MElement& GroupOfDof::getGeoElement(void) const{
  return *element;
}

inline int GroupOfDof::getOrientation(const int dofId) const{
  return (*direction)[dofId];
}

inline void GroupOfDof::orientation(const std::vector<int>& orientation){
  direction = &orientation;
}

#endif
