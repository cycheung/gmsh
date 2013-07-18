#ifndef _DOF_H_
#define _DOF_H_

#include <string>

/**
   @class Dof
   @brief This class represents Degrees Of Freedom

   This class represents Degrees Of Freedom.@n

   A Dof is defined by a pair of to integers called (@c entity, @c type).@n
   By themselfs, these integers have no meaning, they just @em define a Dof.
*/


class Dof{
 private:
  size_t entity;
  size_t type;

 public:
   Dof(void);
   Dof(const Dof& other);
   Dof(size_t entity, size_t type);
  ~Dof(void);

  size_t getEntity(void) const;
  size_t getType(void)   const;

  void setEntity(size_t entity);
  void setType(size_t type);
  void setDof(size_t entity, size_t type);

  bool operator<(const Dof& other) const;
  bool operator>(const Dof& other) const;
  bool operator==(const Dof& other) const;

  std::string toString(void) const;
};


/**
   @fn Dof::Dof(void)
   Instanciates a new Dof with:
   @li @c entity = 0
   @li @c type = 0
   **

   @fn Dof::Dof(const Dof&)
   @param other An other Dof

   Instanciates a @em copy of the given Dof
   **

   @fn Dof::Dof(size_t, size_t)
   @param entity A natural number
   @param type A natural number

   Instanciates a new Dof with the given
   pair (@c entity, @c type)
   **

   @fn Dof::~Dof
   Deletes this Dof
   **

   @fn Dof::getEntity
   @return Returns the associated @a entity of the Dof
   **

   @fn Dof::getType
   @return Returns the associated @a type of the Dof
   **

   @fn Dof::setEntity
   @param entity A natural number

   Sets this Dof @c entity to the given value
   **

   @fn Dof::setType
   @param type A natural number

   Sets this Dof @c type to the given value
   **

   @fn Dof::setDof
   @param entity A natural number
   @param type A natural number

   Sets this Dof to the given pair (@c entity, @c type)
   **

   @fn bool Dof::operator<(const Dof& other) const
   @param other An other Dof to compare the current one
   @return Returns :
   @li @c true if the current Dof is @em smaller than @c other
   @li @c false otherwise
   **

   @fn bool Dof::operator>(const Dof& other) const
   @param other An other Dof to compare the current one
   @return Returns :
   @li @c true if the current Dof is @em greater than @c other
   @li @c false otherwise
   **

   @fn bool Dof::operator==(const Dof& other) const
   @param other An other Dof to compare the current one
   @return Returns :
   @li @c true if the current Dof is @em equal to @c other
   @li @c false otherwise
   **

   @fn Dof::toString
   @return Returns the Dof's string
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline size_t Dof::getEntity(void) const{
  return entity;
}

inline size_t Dof::getType(void) const{
  return type;
}

inline void Dof::setEntity(size_t entity){
  this->entity = entity;
}

inline void Dof::setType(size_t type){
  this->type = type;
}

inline void Dof::setDof(size_t entity, size_t type){
  this->entity = entity;
  this->type   = type;
}

inline bool Dof::operator<(const Dof& other) const{
  return (entity < other.entity) ||
    ((entity == other.entity) && (type < other.type));
}

inline bool Dof::operator>(const Dof& other) const{
  return (entity > other.entity) ||
    ((entity == other.entity) && (type > other.type));
}

inline bool Dof::operator==(const Dof& other) const{
  return (entity == other.entity) && (type == other.type);
}

#endif
