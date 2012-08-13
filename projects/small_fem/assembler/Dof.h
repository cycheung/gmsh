#ifndef _DOF_H_
#define _DOF_H_

#include <string>

/**
   @class Dof
   @brief This class represents degrees of freedom (Dof)
   
   This class represents degrees of freedom (Dof).@n
   These are the terms that will be assembled in the system.

   A Dof is defined by a pair of to integers called (@c entity, @c type).@n
   By themselfs, these integers have no meaning,
   they just @em define a Dof.
*/


class Dof{
 private:
  unsigned int entity;
  unsigned int type;

 public:
   Dof(void);
   Dof(const Dof& other); 
   Dof(unsigned int entity, unsigned int type);
  ~Dof(void);

  unsigned int getEntity(void) const;
  unsigned int getType(void) const;

  void setEntity(unsigned int entity);
  void setType(unsigned int type);
  void setDof(unsigned int entity, unsigned int type);

  bool operator<(const Dof& other) const;
  bool operator>(const Dof& other) const;
  bool operator==(const Dof& other) const;

  std::string toString(void) const;
};


/**
   @fn int Dof::getEntity(void) const
   @return Returns the associated @a entity of the Dof

   @fn int Dof::getType(void) const
   @return Returns the associated @a type of the Dof

   @fn bool Dof::operator<(const Dof& other) const
   @param other An other Dof to compare the current one
   @return Returns :
   @li @c true if the current Dof is @em smaller than @c other
   @li @c false otherwise

   @fn bool Dof::operator>(const Dof& other) const
   @param other An other Dof to compare the current one
   @return Returns :
   @li @c true if the current Dof is @em greater than @c other
   @li @c false otherwise
  
   @fn bool Dof::operator==(const Dof& other) const
   @param other An other Dof to compare the current one
   @return Returns :
   @li @c true if the current Dof is @em equal to @c other
   @li @c false otherwise

   @fn std::string Dof::toString(void) const
   @return Returns the Dof's string
*/

//////////////////////
// Inline Functions //
//////////////////////

inline unsigned int Dof::getEntity(void) const{
  return entity;
}

inline unsigned int Dof::getType(void) const{
  return type;
}

inline void Dof::setEntity(unsigned int entity){
  this->entity = entity;
}

inline void Dof::setType(unsigned int type){
  this->type = type;
}

inline void Dof::setDof(unsigned int entity, unsigned int type){
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
