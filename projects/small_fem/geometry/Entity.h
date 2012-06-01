#ifndef _ENTITY_H_
#define _ENTITY_H_

#include <string>

/**
   @class Entity
   @brief Common methods to all Mesh entities (Node, Edge, ...)
   
   This class represents the common methods to all Mesh entities.
   An entity may be:
   @li A Node
   @li An Edge
   @li A Face
   @li ...

   @warning
   All @em specific mesh entity @em must @em specialize this class. 

   A Entity is defined by multiple parameters:
   @li A Mesh @c ID
   @li A type (Node, Edge, Face, ...)
   @li A value
   @li A physical @c ID in the Mesh

   @note
   Note that the users can't @em instantiated an Entity,
   it is the Mesh job
*/

class Mesh;

class Entity{
 protected:
  int id;
  int type;

  bool   hasValue;
  double value;

  bool   hasPhysical;
  int    physical;

  friend class Mesh;

 protected:
  Entity(const int id, const int type);
  virtual ~Entity(void);

 public:
  bool   gotValue(void) const;
  void   setValue(const double x);
  double getValue(void) const;

  bool   gotPhysical(void) const;
  void   setPhysical(const int physical);
  double getPhysical(void) const;  

  int getId(void) const;
  int getType(void) const;

  bool operator<(const Entity& other) const;
  bool operator==(const Entity& other) const;

  virtual std::string toString(void) const;
};


/**
   @fn Entity::gotValue
   @return Returns:
   @li @c true, if this Entity has a Value assigned
   @li @c false, otherwise

   @fn Entity::setValue
   @param x The value to which we want to set the Entity
   @return Sets the Entity to the given value

   @fn Entity::getValue
   @return Returns the value of the Entity 
   @warning If no value has been set, 
   this method has no specific behaviour

   @fn Entity::gotPhysical
   @return Returns:
   @li @c true, if this Entity is associated with 
   a physical in the Mesh
   @li @c false, otherwise

   @fn Entity::setPhysical
   @param physical The physical @c ID to affect to the Entity
   @return Affect the given physical @c ID to the Entity

   @fn Entity::getPhysical
   @return Returns the physical @c ID associated to the Entity 
   @warning If no physical @c ID has been set, 
   this method has no specific behaviour

   @fn Entity::getId(void)
   @return Returns the @em @c ID of the Entity 
   in the Mesh
  
   @fn Entity::getType(void) const
   @return Returns the @em type of this Entity
   (see description for more details on types)
   
   @fn bool Entity::operator<(const Entity& other) const
   @param other An other Entity
   @return Returns:
   @li @c true, if @em this Entity has a @em smaller
   @em @c ID than the @em other Entity
   @li @c false, otherwise
   @note
   Note that both Entities @em must be of the same type,
   otherwise this operator is allways @em @c flase 

   @fn bool Entity::operator==(const Entity& other) const
   @param other An other Entity
   @return Returns:
   @li @c true, if @em this Entity has the @em same
   @em @c ID as the @em other Entity 
   @li @c false, otherwise
   @note
   Note that both Entities @em must be of the same type,
   otherwise this operator is allways @em @c flase 

   @fn Entity::toString
   @return Returns a string corresponding to the Entity
*/

//////////////////////
// Inline Functions //
//////////////////////

inline Entity::~Entity(void){
}

inline bool Entity::gotValue(void) const{
  return hasValue;
}

inline void Entity::setValue(const double x){
  hasValue = true;
  value = x;
}

inline double Entity::getValue(void) const{
  return value;
}

inline bool Entity::gotPhysical(void) const{
  return hasPhysical;
}

inline void Entity::setPhysical(const int physical){
  hasPhysical = true;
  this->physical = physical;
}

inline double Entity::getPhysical(void) const{
  return physical;
}

inline int Entity::getId(void) const{
  return id;
}

inline int Entity::getType(void) const{
  return type;
}

inline bool Entity::operator<(const Entity& other) const{
  return ((*this).type == other.type) &&
         ((*this).id    < other.id);
}

inline bool Entity::operator==(const Entity& other) const{
  return ((*this).type == other.type) &&
         ((*this).id   == other.id);
}

#endif
