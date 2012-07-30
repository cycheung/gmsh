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
*/

/**
   @class DofComparator
   @brief A class to compare two Dof%s

   A class to compare two Dof%s.
*/

class DofManager;

class Dof{
 private:
  unsigned int entity;
  unsigned int type;

 public:
   Dof(const unsigned int entity, const unsigned int type);
  ~Dof(void);

  unsigned int getEntity(void) const;
  unsigned int getType(void) const;

  bool operator<(const Dof& other) const;
  bool operator>(const Dof& other) const;
  bool operator==(const Dof& other) const;

  std::string toString(void) const;
};

class DofComparator{
 public:
  bool operator()(const Dof* a, const Dof* b) const;
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

   @fn bool DofComparator::operator()(const Dof* a, const Dof* b) const
   @param a A Dof
   @param b Another Dof
   @return operator() is:
   @li @c true, if a is @em smaller than b  
   @li @c false, otherwise
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

inline bool DofComparator::operator()(const Dof* a, const Dof* b) const{
  return *a < *b;
}

#endif
