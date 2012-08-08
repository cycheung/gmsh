#ifndef _GROUPOFELEMENT_H_
#define _GROUPOFELEMENT_H_

#include <string>
#include <vector>

#include "Group.h"
#include "GEntity.h"
#include "Mesh.h"
#include "MElement.h"

#include "GroupOfVertex.h"
#include "GroupOfEdge.h"

/**
   @class GroupOfElement
   @brief A Group of MElement%s

   This class is collection of @em discrete elements (MElement%s).@n
   This class is @em Group.
*/

class Mesh;
class GroupOfVertex;
class GroupOfEdge;

class GroupOfElement: public Group{
 private:
  Mesh*    mesh;
  GEntity* entity;

  static unsigned int nextId;
  unsigned int        id;

  unsigned int            nElement;
  std::vector<MElement*>*  element;

  GroupOfVertex*   gov;
  GroupOfEdge*     goe;

 public:
  GroupOfElement(GEntity& entity, Mesh& mesh);
  virtual ~GroupOfElement(void);

  virtual unsigned int getNumber(void) const;
  virtual unsigned int getId(void) const;
  virtual unsigned int getType(void)   const;

  MElement&                     get(unsigned int i) const;  
  const std::vector<MElement*>& getAll(void) const;  

  GEntity& getEntity(void) const;
  Mesh&    getMesh(void) const;

  GroupOfVertex& getGroupOfVertex(void);
  GroupOfEdge&   getGroupOfEdge(void);

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

inline unsigned int GroupOfElement::getNumber(void) const{
  return nElement;
}

inline unsigned int GroupOfElement::getId(void) const{
  return id;
}

inline unsigned int GroupOfElement::getType(void) const{
  return 1;
}

inline MElement& GroupOfElement::get(unsigned int i) const{
  return *((*element)[i]);
}

inline const std::vector<MElement*>& 
GroupOfElement::getAll(void) const{
  return *element;
}

inline GEntity& GroupOfElement::getEntity(void) const{
  return *entity;
}

inline Mesh& GroupOfElement::getMesh(void) const{
  return *mesh;
}

#endif
