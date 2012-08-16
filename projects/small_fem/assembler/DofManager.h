#ifndef _DOFMANAGER_H_
#define _DOFMANAGER_H_

#include <string>
#include <set>
#include <map>
#include <vector>

#include "FunctionSpace.h"

#include "Dof.h"
#include "GroupOfDof.h"
#include "Comparators.h"

/**
   @class DofManager
   @brief This class manages the Degrees Of Freedom (Dof)
   
   This class manages the Degrees Of Freedom (Dof).@n

   It can map a Dof a @em global @c ID.@n
   Those @c IDs are handeld by the DofManager itself.@n

   If a @em group of Dof got a special meaning 
   (@e e.g: members of the same element),@n
   a DofManager can instantiate (and store) a GroupOfDof.

   Finaly, this class allows to fix a Dof to a given value 
   (in that case, we have a fixed Dof).@n
   Note that we call an @em unknown, a Dof that is @em not fixed.
   

   @warning
   Up to know, a mapped Dof @em can't be @em deleted.@n
   It is also @em impossible to create a @em non @em Mesh @em Related Dof.
*/

class FunctionSpace;

class DofManager{
 private:
  std::set<const Dof*, DofComparator>*         dof;
  std::vector<GroupOfDof*>*                    group;

  std::map<const Dof*, int, DofComparator>*    globalId;
  std::map<const Dof*, double, DofComparator>* fixedDof;

  int nextId;
  
 public:
   DofManager(const FunctionSpace& fs);
  ~DofManager(void);

  int dofNumber(void) const;
  int groupNumber(void) const;

  const std::vector<const Dof*>   getAllDofs(void) const;
  const std::vector<GroupOfDof*>& getAllGroups(void) const;

  int getGlobalId(const Dof& dof) const;

  bool isUnknown(const Dof& dof) const;
  bool fixValue(const Dof& dof, double value);
  std::pair<bool, double> getValue(const Dof& dof) const;

  std::string toString(void) const;

 private:
  void insertDof(Dof& d, GroupOfDof* god);  
};


/**
   @fn DofManager::DofManager
   @param fs A Function Space
   
   Instantiates a new DofManager, based on the given 
   Function Space@n

   This Function Space is used to access:
   @li The Geometric support of the Function Space
   @li The Dof%s associated to each basis function 
   of the Function Space
   **

   @fn DofManager::~DofManager
   @return Deletes this DofManager
   **

   @fn DofManager::dofNumber
   @return Returns the number of Dof in the DofManager
   **

   @fn int DofManager::groupNumber
   @return Returns the number of GroupOfDof in the DofManager   
   **

   @fn DofManager::getAllDofs
   @return Returns all the Dof%s in the DofManager 
   **

   @fn  DofManager::getAllGroups
   @return Returns all the GroupOfDof%s in the DofManager
   **

   @fn DofManager::getGlobalId
   @param dof The Dof from which we want the @em global @c ID
   @return Returns the @em global @em @c ID of the given Dof
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
   @return Returns an std::pair, where:
   <ul>
     <li> The first value is:
       <ol>
         <li> @c true, if the given Dof @em is a @em fixed Dof
         <li> @c false, otherwise
       </ol>
   
     <li> The second value is:
       <ol>
         <li> If the first value was @em @c true, 
	 equal the value of the given (fixed) Dof
         <li> Not specified otherwise
       </ol>
   </ul>
   **

   @fn  DofManager::toString
   @return Returns the DofManager's string
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline int DofManager::dofNumber(void) const{
  return dof->size();
}

inline int DofManager::groupNumber(void) const{
  return group->size();
}

inline const std::vector<const Dof*> DofManager::getAllDofs(void) const{
  return std::vector<const Dof*>(dof->begin(), dof->end());
}

inline const std::vector<GroupOfDof*>& DofManager::getAllGroups(void) const{
  return *group;
}

inline bool DofManager::isUnknown(const Dof& dof) const{
  return fixedDof->count(&dof) == 0;
}

#endif
