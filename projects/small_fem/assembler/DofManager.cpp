#include <sstream>
#include "Exception.h"
#include "GroupOfDof.h"
#include "DofManager.h"

using namespace std;

const unsigned int DofManager::isFixed = 0 - 1; // Largest unsigned int

DofManager::DofManager(void){
  globalIdV = NULL;
  globalIdM = new map<const Dof*, unsigned int, DofComparator>;
  fixedDof  = new map<const Dof*, double, DofComparator>;
}

DofManager::~DofManager(void){
  if(globalIdV){
    for(unsigned int i = 0; i < sizeV; i++)
      delete[] globalIdV[i].id;

    delete[] globalIdV;
  }

  if(globalIdM)
    delete globalIdM;

  delete fixedDof;
}

void DofManager::addToDofManager(const vector<GroupOfDof*>& god){
  // Check if map is still their //
  if(!globalIdM)
    throw
      Exception
      ("DofManager: global id space generated -> can't add Dof");

  // Number Dof //
  const unsigned int nGoD = god.size();

  // Add to DofManager //
  for(unsigned int i = 0; i < nGoD; i++){
    // Dof from god[i]
    const vector<const Dof*>& dof = god[i]->getAll();

    // Init map entry
    const unsigned int nDof = dof.size();

    for(unsigned int j = 0; j < nDof; j++)
      globalIdM->insert(pair<const Dof*, unsigned int>(dof[j], 0));
  }
}

void DofManager::generateGlobalIdSpace(void){
  const map<const Dof*, unsigned int, DofComparator>::iterator
    end = globalIdM->end();

  map<const Dof*, unsigned int, DofComparator>::iterator
    it = globalIdM->begin();

  unsigned int id = 0;

  for(; it != end; it++){
    // Check if unknown
    if(it->second != isFixed){
      it->second = id;
      id++;
    }
  }

  serialize();
  delete globalIdM;
  globalIdM = NULL;
}

void DofManager::serialize(void){
  // Get Data //
  map<const Dof*, unsigned int, DofComparator>::iterator end =
    globalIdM->end();

  map<const Dof*, unsigned int, DofComparator>::iterator it =
    globalIdM->begin();

  // Take the last element *IN* map
  end--;

  first   = it->first->getEntity();
  last    = end->first->getEntity();
  sizeV   = last - first + 1;
  nTotDof = globalIdM->size();

  // Reset 'end': take the first element *OUTSIDE* map
  end++;


  // Alloc //
  globalIdV = new Data[sizeV];


  // Populate //
  map<const Dof*, unsigned int, DofComparator>::iterator
    currentEntity = it;

  // Iterate on vector
  for(unsigned int i = 0; i < sizeV; i++){
    // No dof for this entity
    globalIdV[i].nDof = 0;

    // 'currentEntity - first' matches 'i' ?
    if(it != end && currentEntity->first->getEntity() - first == i)
      // Count all elements with same entity in map
      for(; it !=end &&
            currentEntity->first->getEntity() == it->first->getEntity(); it++)
        globalIdV[i].nDof++; // New Dof found

    // Alloc
    if(globalIdV[i].nDof)
      globalIdV[i].id = new unsigned int[globalIdV[i].nDof];

    else
      globalIdV[i].id = NULL;

    // Add globalIds in vector for this entity
    it = currentEntity;
    for(unsigned int j = 0; j < globalIdV[i].nDof; j++, it++)
      globalIdV[i].id[j] = it->second; // Copy globalId from map

    // Current entity is added to vector:
    //                 go to next entity
    currentEntity = it;
  }
}

pair<bool, unsigned int> DofManager::findSafe(const Dof& dof) const{
  // Is 'dof' in globalIdV range ?
  unsigned int tmpEntity = dof.getEntity();

  if(tmpEntity < first || tmpEntity > last)
    return pair<bool, unsigned int>(false, 42);

  // Offset Entity & Get Type
  const unsigned int entity = tmpEntity - first;
  const unsigned int type   = dof.getType();

  // Look for Entity in globalIdV
  if(globalIdV[entity].nDof > 0 && type <= globalIdV[entity].nDof)
    // If we have Dofs associated to this Entity,
    // get the requested Type and return Id
    return pair<bool, unsigned int>(true, globalIdV[entity].id[type]);

  else
    // If no Dof, return false
    return pair<bool, unsigned int>(false, 42);
}

unsigned int DofManager::getGlobalIdSafe(const Dof& dof) const{
  const pair<bool, unsigned int> search = findSafe(dof);

  if(!search.first)
    throw
      Exception("Their is no Dof %s", dof.toString().c_str());

  else
    return search.second;
}

bool DofManager::isUnknown(const Dof& dof) const{
  if(globalIdM)
    return isUnknownFromMap(dof);

  else
    return isUnknownFromVec(dof);
}

bool DofManager::isUnknownFromVec(const Dof& dof) const{
  const pair<bool, unsigned int> search = findSafe(dof);

  if(search.first)
    return search.second == isFixed;

  else
    throw Exception("Unknown Dof: %s", dof.toString().c_str());
}

bool DofManager::isUnknownFromMap(const Dof& dof) const{
  const map<const Dof*, unsigned int, DofComparator>::
    iterator it = globalIdM->find(&dof);

  if(it != globalIdM->end())
    return it->second == isFixed;

  else
    throw Exception("Unknown Dof: %s", dof.toString().c_str());
}

bool DofManager::fixValue(const Dof& dof, double value){
  // Check if map is still their
  if(!globalIdM)
    throw
      Exception
      ("DofManager: global id space generated -> can't fix values");

  // Get *REAL* Dof
  const map<const Dof*, unsigned int, DofComparator>::iterator it =
    globalIdM->find(&dof);

  // Check if 'dof' exists
  if(it == globalIdM->end())
    return false; // 'dof' doesn't exist

  // If 'dof' exists: it becomes fixed
  fixedDof->insert(pair<const Dof*, double>(it->first, value));
  it->second = isFixed;
  return true;
}

double DofManager::getValue(const Dof& dof) const{
  const map<const Dof*, double, DofComparator>::iterator end =
    fixedDof->end();

  map<const Dof*, double, DofComparator>::iterator it =
    fixedDof->find(&dof);

  if(it == end)
    throw Exception("Dof %s not fixed", dof.toString().c_str());

  return it->second;
}

string DofManager::toString(void) const{
  if(globalIdM)
    return toStringFromMap();

  else
    return toStringFromVec();
}

string DofManager::toStringFromMap(void) const{
  stringstream s;

  const map<const Dof*, unsigned int, DofComparator>::iterator end =
    globalIdM->end();

  map<const Dof*, unsigned int, DofComparator>::iterator it =
    globalIdM->begin();

  for(; it != end; it++){
    s << "("  << it->first->toString() << ": ";

    if(it->second == isFixed)
      s << fixedDof->find(it->first)->second << " -- Fixed value";

    else
      s << it->second                        << " -- Global ID";

    s << ")"  << endl;
  }

  return s.str();
}

string DofManager::toStringFromVec(void) const{
  stringstream s;
  unsigned int nDof;
  pair<bool, unsigned int> search;

  for(unsigned int entity = 0; entity < sizeV; entity++){
    nDof = globalIdV[entity].nDof;

    for(unsigned int type = 0; type < nDof; type++){
      Dof dof(entity + first, type);
      search = findSafe(dof);

      if(search.first){
        s << "("  << dof.toString() << ": ";

        if(search.second == isFixed)
          s << fixedDof->find(&dof)->second << " -- Fixed value";

        else
          s << search.second                << " -- Global ID";

        s << ")"  << endl;
      }
    }
  }

  return s.str();
}
