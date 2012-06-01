#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include <string>
#include <vector>
#include "Entity.h"
#include "Node.h"
#include "Jacobian.h"

/**
   @class Element
   @brief Represents a geometrical element of a Mesh

   An Element is a geometrical element of a Mesh.@n

   An Element is defined by:
   @li An @c ID in the Mesh
   @li A set of @em Node%s, giving the @em position of the 
   Element in the Mesh
   @li A set of @em Entity, giving the @em geomtrical
   @em elements of the Element
   @li A set of @em Orientations for the Entity
   @li A type
   
   @note
   In the case of @em nodal elements, the sets of
   @em Node%s and @em Entity are @em equal

   Element%s type are given by the @em concatenation of:
   @li The number of @em Node%s describing the Element

   And

   @li @c 0 for @em nodal Element%s
   @li @c 1 for @em edge Element%s

   @warning This may change in the futur
   @todo Find a clever way for Element%s types

   Moreover, an Element will instantiate a Jacobian 
   transformation.

   @note
   An Element can't be instantiated by the users,
   it is the @em Mesh job
 */

class Mesh;

class Element{
 private:
  int type;
  int id;
  int numberOfEntity;

  std::vector<Node*>* node;

  std::vector<Entity*>* entity;
  std::vector<int>* orientation;
  std::vector<int>* entityId;  

  Jacobian* jac;

  friend class Mesh;

 private:
   Element(const int id, const int type);
  ~Element(void);

 public:
  Entity& getEntity(const int id) const;
  const std::vector<Entity*>& getAllEntities(void) const;
  const std::vector<int>&   getAllEntitiesId(void);
  
  const std::vector<Node*>& getAllNodes(void) const;

  int getOrientation(const int id) const;
  const std::vector<int>& getAllOrientations(void) const;  

  int getId(void) const;
  int nEntity(void) const;
  int getType(void) const;
  
  Jacobian& getJacobian(void) const;

  bool operator==(const Element& other) const;
  bool operator!=(const Element& other) const;
    
  std::string toString(void) const;

 private:
  void buildJacobian(void);
};


/**
   @fn Element::getEntity
   @param id The Entity @em local @c ID in the Element
   @return Returns the requested Entity
   @note
   Entity @em local @c IDs are ranging from 
   @em 0 to @em E @em -  @em 1,@n
   where E is the number of Entity in the Element

   @fn Element::getAllEntities
   @return Returns the set of all the Entity of the Element

   @fn Element::getAllEntitiesId
   @return Returns the set of all the Entity @c IDs 
   of the Element
   @warning
   The @c IDs are related to the Mesh, so they are
   @em global @c IDs

   @fn Element::getAllNodes
   @return Returns the set of all the Node%s of the Element 

   @fn Element::getOrientation
   @param id An Element Entity @em local @c ID
   @return Returns the orientation of the given 
   (by its @em local @c ID) Entity

   @fn Element::getAllOrientations
   @return Returns the set of orientations of the  
   Element Entity
   
   @fn Element::getId
   @return Returns the @c ID of this Element in the Mesh
   
   @fn Element::nEntity
   @return Returns the @em size of set of Entity

   @fn Element::getType
   @return Returns the @em type of this Element

   @fn Element::getJacobian
   @return Returns the Jacobian allowing
   the mapping of this Element onto
   a reference element

   @fn bool Element::operator==(const Element& other) const;
   @param other An other Element
   @return Returns:
   @li @c true, if this Element and the other are 
   the @em same
   @li @c false, otherwise
   
   @fn bool Element::operator!=(const Element& other) const;
   @param other
   @return Returns:
   @li @c true, if this Element and the other are 
   @em not the @em same
   @li @c false, otherwise
   
   @fn Element::toString
   @return Returns a string associated to this Element

*/

//////////////////////
// Inline Functions //
//////////////////////

inline Entity& Element::getEntity(const int id) const{
  return *((*entity)[id]);
}

inline const std::vector<Entity*>& Element::getAllEntities(void) const{
  return *entity;
}

inline const std::vector<Node*>& Element::getAllNodes(void) const{
  return *node;
}

inline int Element::getOrientation(const int id) const{
  return (*orientation)[id];
}

inline const std::vector<int>& Element::getAllOrientations(void) const{
  return *orientation;
}

inline int Element::getId(void) const{
  return id;
}

inline int Element::nEntity(void) const{
  return numberOfEntity;
}

inline int Element::getType(void) const{
  return type;
}

inline Jacobian& Element::getJacobian(void) const{
  return *jac;
}

inline bool Element::operator==(const Element& other) const{
  return (*this).id == other.id;
}

inline bool Element::operator!=(const Element& other) const{
  return !(*this == other);
}

#endif
