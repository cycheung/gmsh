#ifndef _DOFMANAGER_H_
#define _DOFMANAGER_H_

#include <string>
#include <set>
#include <map>
#include <vector>

#include "FunctionSpace.h"
#include "MElement.h"
#include "GroupOfDof.h"

/**
   @class DofManager
   @brief This class manages the degrees freedom (Dof)
   
   This class manages the degrees freedom (Dof).@n

   It can map a Dof a @em global @c ID.@n
   Those @c IDs are handeld by the DofManager itself.@n

   It can also map a given Dof to the corresponding @em Entity in the Mesh.@n

   A DofManager can be instantiated from a list of @em Element%s.@n

   If a @em group of Dof got a special meaning 
   (@e e.g: members of the same element),@n
   a DofManager can instantiate (and map) GroupOfDof%s.

   @note
   A DofManager is the @em only @em one allowed to instatiate a Dof.@n
   It is also the @em only @em one allowed to instatiate a GroupOfDof.

   @warning
   Up to know, a Dof @em can't be @em deleted.@n
   It is also @em impossible to create a @em non @em Mesh @em Related Dof.

   @todo
   A more @em general DofManager, with non Mesh Dof, etc@n
   Allow hybrid mesh
*/

/**
   @class DofManager::DofComparator
   @brief A private class of DofManager, that can compare two Dof%s

   A private class of DofManager, that can compare two Dof%s.
*/

class DofManager{
 private:

  class DofComparator{
  public:
    bool operator()(const Dof* a, const Dof* b) const;
  };

  const FunctionSpace* fs;

  std::set<Dof*, DofComparator>*      dof;
  std::vector<GroupOfDof*>*           group;
  std::map<Dof*, int, DofComparator>* globalId;

  int nTotVertex;
  int nextId;
  
 public:
   DofManager(const FunctionSpace& fs);
  ~DofManager(void);

  int dofNumber(void) const;
  int groupNumber(void) const;

  const std::vector<Dof*>         getAllDofs(void) const;
  const std::vector<GroupOfDof*>& getAllGroups(void) const;

  int getGlobalId(Dof& dof) const;

  void setAsConstant(const GroupOfElement& goe, double value);

  std::string toString(void) const;

 private:
  Dof** dofFromElement(MElement& element, int* nDof);
  void  insertDof(Dof* d, GroupOfDof* god);
  
};


/**
   @fn DofManager::DofManager(const std::vector<MElement*>& element)
   @param element A list of Element%s from which Dof%s will be instantiated
   @return Instantiate a new DofManager from the given list of Element%s

   @fn DofManager::~DofManager(void)
   @return Deletes the DofManager

   @fn int DofManager::dofNumber(void) const
   @return Returns the number of Dof mapped in the DofManager

   @fn int DofManager::groupNumber(void) const
   @return Returns the number of GroupOfDof mapped in the DofManager   

   @fn const std::vector<Dof*>& DofManager::getAllDofs(void) const
   @return Returns all the Dof%s in the DofManager 

   @fn const std::multimap<int, Dof*>& DofManager::getAllPhysicals(void) const
   @return Returns a map between @e Physical @e @c IDs and Dof%s
   @warning This should be replaced with a method @em hiding the multimap
   @todo Replace with a method @em hiding the multimap

   @fn const std::vector<GroupOfDof*>& DofManager::getAllGroups(void) const
   @return Returns all the GroupOfDof%s in the DofManager

   @fn int DofManager::getGlobalId(Dof& dof) const
   @param dof The Dof from which we want the @em global @c ID
   @return Returns the @em global @em @c ID of the given Dof
 
   @fn Entity& DofManager::getEntity(Dof& dof) const
   @param dof The Dof from which we want the @em entity
   @return Returns the @em entity associated with the given Dof

   @fn std::string DofManager::toString(void) const
   @return Returns the DofManager's string

   @fn bool DofManager::DofComparator::operator()(const Dof* a, const Dof* b) const
   @param a A Dof
   @param b Another Dof
   @return operator() is:
   @li @c true, if a is @em smaller than b  
   @li @c false, otherwise
}
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

inline const std::vector<Dof*> DofManager::getAllDofs(void) const{
  return std::vector<Dof*>(dof->begin(), dof->end());
}

inline const std::vector<GroupOfDof*>& DofManager::getAllGroups(void) const{
  return *group;
}

inline int DofManager::getGlobalId(Dof& dof) const{
  return globalId->find(&dof)->second;
}

inline bool DofManager::DofComparator::operator()(const Dof* a, const Dof* b) const{
  return *a < *b;
}

#endif
