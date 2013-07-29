#ifndef _DOFMANAGER_H_
#define _DOFMANAGER_H_

#include <string>
#include <vector>
#include <map>

#include "GroupOfDof.h"
#include "Comparators.h"

/**
   @class DofManager
   @brief This class manages the Degrees of Freedom (Dof)

   This class numbers the Degrees of Freedom (Dof).

   It can map a Dof to a unique number, called global ID.

   In addtion, this class allows to assign a Dof to a given value.
   A Dof that has been assigned to a value is called a fixed Dof.
   The global ID of a fixed Dof is not unique and is equal to
   DofManager::isFixedId().

   Finaly, the global IDs given to the unfixed Dof%s ranges from 0 to
   (total number of Dof - number of fixed Dof - 1).
*/

class DofManager{
 private:
  static const size_t isFixed;

  std::vector<std::vector<size_t> > globalIdV;
  std::map<Dof, size_t>             globalIdM;
  std::map<Dof, double>             fixedDof;

  size_t first;
  size_t last;
  size_t nTotDof;

 public:
  static const size_t isFixedId(void);

 public:
   DofManager(void);
  ~DofManager(void);

  void addToDofManager(const std::vector<GroupOfDof*>& god);
  void generateGlobalIdSpace(void);

  size_t getGlobalId(const Dof& dof)     const;
  size_t getGlobalIdSafe(const Dof& dof) const;

  bool   fixValue(const Dof& dof, double value);
  double getValue(const Dof& dof) const;

  size_t getTotalDofNumber(void) const;
  size_t getUnfixedDofNumber(void) const;
  size_t getFixedDofNumber(void) const;

  std::string toString(void) const;

 private:
  void serialize(void);
  std::pair<bool, size_t> findSafe(const Dof& dof) const;

  std::string toStringFromMap(void) const;
  std::string toStringFromVec(void) const;
};


/**
   @fn DofManager::isFixedId

   Fixed Dof got a special global ID (which is the highest possible number).

   @return The special ID of fixed Dof
   **

   @fn DofManager::DofManager

   Instantiates a new DofManager
   **

   @fn DofManager::~DofManager

   Deletes this DofManager
   **

   @fn DofManager::addToDofManager
   @param god A vector of GroupOfDof

   Adds the given Dof%s in this DofManager.
   The same Dof may be insterd multiple time,
   but it will be given the same unique ID.
   **

   @fn DofManager::generateGlobalIdSpace

   Numbers every non fixed Dof of this DofManager.
   Each Dof will be given a unique global ID.
   **

   @fn DofManager::getGlobalId
   @param dof The Dof from which we want the global ID
   @return Returns the global ID of the given Dof

   If the requested Dof has not been added to this DofManager,
   or if DofManager::generateGlobalIdSpace() has not been called,
   the behaviour of this method is unpredicable.
   Actually, it will most probably lead to a process crash.

   @see DofManager::getGlobalIdSafe
   **

   @fn DofManager::getGlobalIdSafe
   @param dof The Dof from which we want the global ID
   @return Returns the global ID of the given Dof

   If the requested Dof has not been added to this DofManager,
   or if DofManager::generateGlobalIdSpace() has not been called,
   an Exception is thrown.

   This method is safer but slower than DofManager::getGlobalId().

   @see DofManager::getGlobalId
   **

   @fn DofManager::fixValue
   @param dof A Dof
   @param value A real number

   Fixes the given Dof to the given value.

   @return Returns:
   @li true, if the operation is a success
   @li false, otherwise

   Here are three important cases, where DofManager::fixValue() will fail:
   @li The given Dof is not in this DofManager
   @li The given Dof is already fixed
   @li DofManager::generateGlobalIdSpace() has been called
   **

   @fn DofManager::getValue
   @param dof A Dof
   @return Returns the value of the given Dof, if it has been fixed

   This method throws an Exception if the Dof has not been fixed
   **

   @fn DofManager::getTotalDofNumber
   @return Returns the number of Dof%s in this DofManager (fixed and unfixed)
   **

   @fn DofManager::getUnfixedDofNumber
   @return Returns the number of Unfixed Dof%s in this DofManager
   **

   @fn DofManager::getFixedDofNumber
   @return Returns the number of fixed Dof%s in this DofManager
   **

   @fn  DofManager::toString
   @return Returns the DofManager string
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline const size_t DofManager::isFixedId(void){
  return isFixed;
}

inline size_t DofManager::getGlobalId(const Dof& dof) const{
  return globalIdV[dof.getEntity() - first][dof.getType()];
}

#endif
