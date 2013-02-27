#include <sstream>
#include "Exception.h"
#include "GroupOfDof.h"
#include "DofFixedException.h"
#include "DofManager.h"

using namespace std;

DofManager::DofManager(void){
  // Alloc Struct //
  globalId     = new map<const Dof*, DofData, DofComparator>;
  nNotFixedDof = 0;
}

DofManager::~DofManager(void){
  delete globalId;
}

void DofManager::addToDofManager(const vector<GroupOfDof*>& god){
  // Number Dof //
  const unsigned int nGoD = god.size();

  // Tmp //
  const DofData tmp = {false, 0, 0};
  std::pair<
    map<const Dof*, DofData, DofComparator>::iterator,
    bool>
    inserted;

  for(unsigned int i = 0; i < nGoD; i++){
    // Dof from god[i]
    const vector<const Dof*>& dof = god[i]->getAll();

    // Init map entry
    const unsigned int nDof = dof.size();

    for(unsigned int j = 0; j < nDof; j++){
      inserted = globalId->insert(pair<const Dof*, DofData>(dof[j], tmp));

      if(inserted.second)
        nNotFixedDof++;
    }
  }
}
/*
void DofManager::generateGlobalIdSpace(bool withFixedValue){
  if(withFixedValue)
    generateGlobalIdSpace();

  else
    generateGlobalIdSpaceWithoutFixedDof();
}

void DofManager::generateGlobalIdSpace(void){
  const map<const Dof*, DofData, DofComparator>::iterator
    end = globalId->end();

  map<const Dof*, DofData, DofComparator>::iterator
    it = globalId->begin();

  for(unsigned int id = 0; it != end; it++, id++)
    it->second.globalId = id;
}

void DofManager::generateGlobalIdSpaceWithoutFixedDof(void){
  const map<const Dof*, DofData, DofComparator>::iterator
    end = globalId->end();

  map<const Dof*, DofData, DofComparator>::iterator
    it = globalId->begin();

  unsigned int id = 0;

  for(; it != end; it++){
    // Check if unknown
    if(!it->second.isFixed){
      it->second.globalId = id;
      id++;
    }
  }
}
*/

void DofManager::generateGlobalIdSpace(void){
  const map<const Dof*, DofData, DofComparator>::iterator
    end = globalId->end();

  map<const Dof*, DofData, DofComparator>::iterator
    it = globalId->begin();

  unsigned int id = 0;

  for(; it != end; it++){
    // Check if unknown
    if(!it->second.isFixed){
      it->second.globalId = id;
      id++;
    }
  }
}

unsigned int DofManager::getGlobalId(const Dof& dof) const{
  const map<const Dof*, DofData, DofComparator>::
    iterator it = globalId->find(&dof);

  if(it == globalId->end())
    throw
      Exception("Their is no Dof %s", dof.toString().c_str());

  else if(it->second.isFixed)
    throw
      DofFixedException(*it->first, it->second.fixedValue);

  else
    return it->second.globalId;
}

bool DofManager::fixValue(const Dof& dof, double value){
  // Get *REAL* Dof
  map<const Dof*, DofData, DofComparator>::iterator it =
    globalId->find(&dof);

  // Check if 'dof' exists
  if(it == globalId->end())
    return false; // 'dof' doesn't exist

  // Remove one Dof from nDofFixedDof
  //  If it wasn't a fixed one !
  if(!it->second.isFixed)
    nNotFixedDof--;

  // Dof becomes fixed
  it->second.isFixed    = true;
  it->second.fixedValue = value;

  return true;
}

pair<bool, double> DofManager::getValue(const Dof& dof) const{
  map<const Dof*, DofData, DofComparator>::iterator end =
    globalId->end();

  map<const Dof*, DofData, DofComparator>::iterator it =
    globalId->find(&dof);

  if(it == end || !it->second.isFixed)
    return pair<bool, double>(false, 42);

  else
    return pair<bool, double>(true, it->second.fixedValue);
}

string DofManager::toString(void) const{
  stringstream s;

  const map<const Dof*, DofData, DofComparator>::iterator end =
    globalId->end();

  map<const Dof*, DofData, DofComparator>::iterator it =
    globalId->begin();

  for(; it != end; it++){
    s << "("  << it->first->toString() << ": ";

    if(it->second.isFixed)
      s << it->second.fixedValue << " -- Fixed value";
    else
      s << it->second.globalId   << " -- Global ID";

    s << ")"  << endl;
  }

  return s.str();
}
