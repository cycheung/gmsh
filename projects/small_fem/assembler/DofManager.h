#ifndef _DOFMANAGER_H_
#define _DOFMANAGER_H_

#include <string>
#include <map>

#include "Dof.h"
#include "Comparators.h"
#include "FunctionSpace.h"

/**
   @class DofManager
   @brief This class manages the Degrees Of Freedom (Dof)

   This class numbers the Degrees Of Freedom (Dof).@n

   It can map a Dof to a @em global @c ID.@n
   Those @c IDs are handeld by the DofManager itself.@n

   Finaly, this class allows to fix a Dof to a given value
   (in that case, we have a fixed Dof).@n
   Note that we call an @em unknown, a Dof that is @em not fixed.

   @warning
   Up to know, a mapped Dof @em can't be @em deleted.@n
*/

class DofManager{
  static const unsigned int isFixed;

  std::map<const Dof*, unsigned int, DofComparator>* globalId;
  std::map<const Dof*, double, DofComparator>*       fixedDof;

 public:
  static const unsigned int isFixedId(void);

 public:
   DofManager(void);
  ~DofManager(void);

  void addToDofManager(const std::vector<GroupOfDof*>& god);
  void generateGlobalIdSpace(void);

  unsigned int getGlobalId(const Dof& dof) const;

  bool   isUnknown(const Dof& dof) const;
  bool   fixValue(const Dof& dof, double value);
  double getValue(const Dof& dof) const;

  unsigned int getDofNumber(void) const;

  std::string toString(void) const;

 private:
  void serialize(void);
};


/**
   @fn DofManager::isFixedId

   Fixed Dof got a special global @c ID.@n
   This global @c ID is returned by this class method.

   @see DofManager::fixValue()
   @see DofManager::getGlobalId()

   @return The special @c ID for fixed Dof
   **

   @fn DofManager::DofManager

   Instantiates a new DofManager
   **

   @fn DofManager::~DofManager

   Deletes this DofManager
   **

   @fn DofManager::addToGlobalIdSpace
   @param god A vector of GroupOfDof

   Adds the given Dof%s in this DofManager
   **

   @fn DofManager::generateGlobalIdSpace

   Numbers every non fixed Dof%s of this DofManager.@n
   Each Dof%s will be given a @em unique @em global @c ID .
   **

   @fn DofManager::getGlobalId
   @param dof The Dof from which we want the @em global @c ID
   @return Returns the @em global @em @c ID of the given Dof

   @note
   If the given Dof has been fixed,
   DofManager::isFixedId is returned

   @note
   If the given Dof is not in this DofManager,
   an Exception is thrown
   **

   @fn DofManager::isUnknown
   @param dof A Dof
   @return Returns:
   @li @c true, if the given Dof is an unknwon
   (@em i.e. a non fixed Dof)
   @li @c false, otherwise
   **

   @fn DofManager::fixValue
   @param dof A Dof
   @param value A real number

   Fixes the given Dof to the given value

   @return Returns:
   @li @c true, if the operation is a success
   @li @c false, otherwise

   @note
   Here are two important cases, where fixValue() will fail:
   @li The given Dof is not in this DofManager
   @li The given Dof is already fixed
   **

   @fn DofManager::getValue
   @param dof A Dof
   @return
   <ul>
     <li> Returns the value of the given Dof, if it has been fixed
     <li> Throws an Exception if the Dof has not been fixed
   </ul>
   **

   @fn DofManager::getDofNumber
   @return Returns the number of Dof%s in
   this DofManager (without fixed Dof%s)
   **

   @fn  DofManager::toString
   @return Returns the DofManager's string
*/

//////////////////////
// Inline Functions //
//////////////////////

inline const unsigned int DofManager::isFixedId(void){
  return isFixed;
}

inline bool DofManager::isUnknown(const Dof& dof) const{
  return globalId->find(&dof)->second == isFixed;
}

inline unsigned int DofManager::getDofNumber(void) const{
  return globalId->size() - fixedDof->size();
}

#endif
