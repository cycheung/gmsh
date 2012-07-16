#ifndef _GROUPOFELEMENTS_H_
#define _GROUPOFELEMENTS_H_

#include <string>
#include <vector>

#include "Group.h"
#include "GEntity.h"
#include "MElement.h"

/**
   @class GroupOfElements
   @brief A Group of MElement%s

   This class is collection of @em discrete elements (MElement%s).@n
   This class is @em Group.
*/


class GroupOfElements: public Group{
 private:
  int                      id;      
  GEntity*                 entity;

  unsigned int            nElement;
  std::vector<MElement*>*  element;

 public:
  GroupOfElements(GEntity& entity, int id);
  virtual ~GroupOfElements(void);

  virtual int getNbElements(void) const;
  virtual MElement& getElement(int i) const;
  
  virtual const std::vector<MElement*>& 
    getAllElements(void) const;

  int      getId(void) const; //! @todo Put in Common interface (Group) ?
  GEntity& getEntity(void) const;

  virtual std::string toString(void) const;
};


/**
   @fn GroupOfElements::GroupOfElements
   Instantiates a new GroupOfElements 
   based on a given GEntity and with the unique ID '@c id' 
   @param entity The GEntity that describes the Group
   @param id A @em unique number 

   @fn GroupOfElements::~GroupOfElements
   Deletes this GroupOfElements
 
   @fn GroupOfElements::getEntity
   @return Returns the Entity used to build 
   this GroupOfElements
*/


//////////////////////
// Inline Functions //
//////////////////////

inline int GroupOfElements::getNbElements(void) const{
  return nElement;
}

inline MElement& GroupOfElements::getElement(int i) const{
  return *((*element)[i]);
}

inline const std::vector<MElement*>& 
GroupOfElements::getAllElements(void) const{
  return *element;
}

inline GEntity& GroupOfElements::getEntity(void) const{
  return *entity;
}

inline int GroupOfElements::getId(void) const{
  return id;
}

#endif
