#ifndef _GROUPOFELEMENTS_H_
#define _GROUPOFELEMENTS_H_

#include <string>
#include <vector>

#include "GroupTyped.h"
#include "GEntity.h"
#include "MElement.h"

/**
   @class GroupOfElements
   @brief A Group of MElement%s

   This class is collection of @em discrete elements (MElement%s).@n
   This class is @em GroupTyped<MElement>.
*/


class GroupOfElements: public GroupTyped<MElement>{
 private:
  int                      id;      
  GEntity*                 entity;

  unsigned int            nElement;
  std::vector<MElement*>*  element;

 public:
  GroupOfElements(GEntity& entity, int id);
  virtual ~GroupOfElements(void);

  virtual int getNumber(void) const;
  virtual int getId(void) const;

  virtual MElement&                     get(int i) const;  
  virtual const std::vector<MElement*>& getAll(void) const;  

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

inline int GroupOfElements::getNumber(void) const{
  return nElement;
}

inline int GroupOfElements::getId(void) const{
  return id;
}

inline MElement& GroupOfElements::get(int i) const{
  return *((*element)[i]);
}

inline const std::vector<MElement*>& 
GroupOfElements::getAll(void) const{
  return *element;
}

inline GEntity& GroupOfElements::getEntity(void) const{
  return *entity;
}

#endif
