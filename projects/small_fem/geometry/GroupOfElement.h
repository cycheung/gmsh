#ifndef _GROUPOFELEMENT_H_
#define _GROUPOFELEMENT_H_

#include <string>
#include <vector>

#include "Group.h"
#include "GEntity.h"
#include "MElement.h"

/**
   @class GroupOfElement
   @brief A Group of MElement%s

   This class is collection of @em discrete elements (MElement%s).@n
   This class is @em Group.
*/


class GroupOfElement: public Group{
 private:
  int                      id;      
  GEntity*                 entity;

  unsigned int            nElement;
  std::vector<MElement*>*  element;

 public:
  GroupOfElement(GEntity& entity, int id);
  virtual ~GroupOfElement(void);

  virtual int getNumber(void) const;
  virtual int getId(void)     const;
  virtual int getType(void)   const;

  MElement&                     get(int i) const;  
  const std::vector<MElement*>& getAll(void) const;  

  GEntity& getEntity(void) const;

  virtual std::string toString(void) const;
};


/**
   @fn GroupOfElement::GroupOfElement
   Instantiates a new GroupOfElement 
   based on a given GEntity and with the unique ID '@c id' 
   @param entity The GEntity that describes the Group
   @param id A @em unique number 

   @fn GroupOfElement::~GroupOfElement
   Deletes this GroupOfElement

   @fn GroupOfElement::get
   @param i An interger ranging from 0 
   to GroupOfElement::getNumber() - 1
   @return Returns the ith element of the Group

   @fn GroupOfElement::getAll
   @return Returns all the elements of the Group
 
   @fn GroupOfElement::getEntity
   @return Returns the Entity used to build 
   this GroupOfElement
*/


//////////////////////
// Inline Functions //
//////////////////////

inline int GroupOfElement::getNumber(void) const{
  return nElement;
}

inline int GroupOfElement::getId(void) const{
  return id;
}

inline int GroupOfElement::getType(void) const{
  return 1;
}

inline MElement& GroupOfElement::get(int i) const{
  return *((*element)[i]);
}

inline const std::vector<MElement*>& 
GroupOfElement::getAll(void) const{
  return *element;
}

inline GEntity& GroupOfElement::getEntity(void) const{
  return *entity;
}

#endif
