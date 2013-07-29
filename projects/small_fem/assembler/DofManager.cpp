#include <sstream>
#include "Exception.h"
#include "GroupOfDof.h"
#include "DofManager.h"

using namespace std;

const size_t DofManager::isFixed = 0 - 1; // Largest size_t

DofManager::DofManager(void){
}

DofManager::~DofManager(void){
}

void DofManager::addToDofManager(const vector<GroupOfDof*>& god){
  // Check if map is still their //
  if(!globalIdM.empty())
    throw
      Exception
      ("DofManager: global id space generated -> can't add Dof");

  // Number Dof //
  const size_t nGoD = god.size();

  // Add to DofManager //
  for(size_t i = 0; i < nGoD; i++){
    // Dof from god[i]
    const vector<Dof>& dof = god[i]->getDof();

    // Init map entry
    const size_t nDof = dof.size();

    for(size_t j = 0; j < nDof; j++)
      globalIdM.insert(pair<Dof, size_t>(dof[j], 0));
  }
}

void DofManager::generateGlobalIdSpace(void){
  map<Dof, size_t>::iterator end = globalIdM.end();
  map<Dof, size_t>::iterator it  = globalIdM.begin();

  size_t id = 0;

  for(; it != end; it++){
    // Check if unknown
    if(it->second != isFixed){
      it->second = id;
      id++;
    }
  }

  serialize();
  globalIdM.clear();
}

void DofManager::serialize(void){
  // Get Data //
  map<Dof, size_t>::iterator end = globalIdM.end();
  map<Dof, size_t>::iterator it  = globalIdM.begin();

  // Take the last element *IN* map
  end--;

  first   = it->first.getEntity();
  last    = end->first.getEntity();
  nTotDof = globalIdM.size();

  // Reset 'end': take the first element *OUTSIDE* map
  end++;


  // Alloc //
  const size_t sizeV = last - first + 1;
  globalIdV.resize(sizeV);

  // Populate //
  size_t nDof;
  map<Dof, size_t>::iterator currentEntity = it;

  // Iterate on vector
  for(size_t i = 0; i < sizeV; i++){
    // No dof found
    nDof = 0;

    // 'currentEntity - first' matches 'i' ?
    if(it != end && currentEntity->first.getEntity() - first == i)
      // Count all elements with same entity in map
      for(; it !=end &&
            currentEntity->first.getEntity() == it->first.getEntity(); it++)
        nDof++; // New Dof found

    // Alloc
    if(nDof)
      globalIdV[i].resize(nDof);

    // Add globalIds in vector for this entity
    it = currentEntity;
    for(size_t j = 0; j < nDof; j++, it++)
      globalIdV[i][j] = it->second; // Copy globalId from map

    // Current entity is added to vector:
    //                 go to next entity
    currentEntity = it;
  }
}

pair<bool, size_t> DofManager::findSafe(const Dof& dof) const{
  // Is globabId Vector allocated ?-
  if(globalIdV.empty())
    throw
      Exception
      ("Cannot get Dof %s ID, since global ID space has not been generated",
       dof.toString().c_str());

  // Is 'dof' in globalIdV range ?
  size_t tmpEntity = dof.getEntity();

  if(tmpEntity < first || tmpEntity > last)
    return pair<bool, size_t>(false, 42);

  // Offset Entity & Get Type
  const size_t entity = tmpEntity - first;
  const size_t type   = dof.getType();

  // Look for Entity in globalIdV
  const size_t nDof = globalIdV[entity].size();

  if(nDof > 0 && type <= nDof)
    // If we have Dofs associated to this Entity,
    // get the requested Type and return Id
    return pair<bool, size_t>(true, globalIdV[entity][type]);

  else
    // If no Dof, return false
    return pair<bool, size_t>(false, 42);
}

size_t DofManager::getGlobalIdSafe(const Dof& dof) const{
  const pair<bool, size_t> search = findSafe(dof);

  if(!search.first)
    throw
      Exception("Their is no Dof %s", dof.toString().c_str());

  else
    return search.second;
}

bool DofManager::fixValue(const Dof& dof, double value){
  // Check if map is still their
  if(globalIdM.empty())
    return false;

  // Get *REAL* Dof
  const map<Dof, size_t>::iterator it = globalIdM.find(dof);

  // Check if 'dof' exists
  if(it == globalIdM.end())
    return false; // 'dof' doesn't exist

  // If 'dof' exists: it becomes fixed
  fixedDof.insert(pair<Dof, double>(it->first, value));
  it->second = isFixed;
  return true;
}

double DofManager::getValue(const Dof& dof) const{
  const map<Dof, double>::const_iterator end = fixedDof.end();
  const map<Dof, double>::const_iterator  it = fixedDof.find(dof);

  if(it == end)
    throw Exception("Dof %s not fixed", dof.toString().c_str());

  return it->second;
}

size_t DofManager::getTotalDofNumber(void) const{
  if(!globalIdM.empty())
    return globalIdM.size();

  else
    return nTotDof;
}

size_t DofManager::getUnfixedDofNumber(void) const{
  if(!globalIdM.empty())
    return globalIdM.size() - fixedDof.size();

  else
    return nTotDof - fixedDof.size();
}

size_t DofManager::getFixedDofNumber(void) const{
  return fixedDof.size();
}

string DofManager::toString(void) const{
  if(!globalIdM.empty())
    return toStringFromMap();

  else
    return toStringFromVec();
}

string DofManager::toStringFromMap(void) const{
  stringstream s;

  map<Dof, size_t>::const_iterator end = globalIdM.end();
  map<Dof, size_t>::const_iterator it  = globalIdM.begin();

  for(; it != end; it++){
    s << "("  << it->first.toString() << ": ";

    if(it->second == isFixed)
      s << fixedDof.find(it->first)->second << " -- Fixed value";

    else
      s << it->second                       << " -- Global ID";

    s << ")"  << endl;
  }

  return s.str();
}

string DofManager::toStringFromVec(void) const{
  const size_t sizeV = globalIdV.size();

  stringstream s;
  size_t nDof;
  pair<bool, size_t> search;

  for(size_t entity = 0; entity < sizeV; entity++){
    nDof = globalIdV[entity].size();

    for(size_t type = 0; type < nDof; type++){
      Dof dof(entity + first, type);
      search = findSafe(dof);

      if(search.first){
        s << "("  << dof.toString() << ": ";

        if(search.second == isFixed)
          s << fixedDof.find(dof)->second << " -- Fixed value";

        else
          s << search.second              << " -- Global ID";

        s << ")"  << endl;
      }
    }
  }

  return s.str();
}
