#ifndef _DOF_H_
#define _DOF_H_

#include <string>

/**
   @class Dof
   @brief This class represents degrees of freedom (Dof)
   
   This class represents degrees of freedom (Dof).@n
   These are the terms that will be assembled in the system.

   A Dof is defined by a pair of to integers called (@c entity, @c type).@n
   By themselfs, these integers have no meaning.@n
   They just @em define a Dof.

   @note
   Note that users are not allowed to instanciate a Dof.@n
   This is the DofManager's responsability.
*/

class DofManager;

class Dof{
 private:
  bool unknown;

  unsigned int entity;
  unsigned int type;

  double value;

  friend class DofManager;

 private:
  // Construction and destruction are not for the user responsablity
   Dof(const unsigned int entity, const unsigned int type);
  ~Dof(void);

 public:
  bool isUnknown(void) const;

  unsigned int getEntity(void) const;
  unsigned int getType(void) const;

  double getValue(void) const;

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

inline bool Dof::isUnknown(void) const{
  return unknown;
}

inline unsigned int Dof::getEntity(void) const{
  return entity;
}

inline unsigned int Dof::getType(void) const{
  return type;
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
