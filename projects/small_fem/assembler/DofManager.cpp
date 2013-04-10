#include <sstream>
#include "Exception.h"
#include "GroupOfDof.h"
#include "DofManager.h"

using namespace std;

const unsigned int DofManager::isFixed = 0 - 1; // Largest unsigned int

DofManager::DofManager(void){
  globalId = new map<const Dof*, unsigned int, DofComparator>;
  fixedDof = new map<const Dof*, double, DofComparator>;
}

DofManager::~DofManager(void){
  delete globalId;
  delete fixedDof;
}

void DofManager::addToDofManager(const vector<GroupOfDof*>& god){
  // Number Dof //
  const unsigned int nGoD = god.size();

  // Add to DofManager //
  for(unsigned int i = 0; i < nGoD; i++){
    // Dof from god[i]
    const vector<const Dof*>& dof = god[i]->getAll();

    // Init map entry
    const unsigned int nDof = dof.size();

    for(unsigned int j = 0; j < nDof; j++){
      globalId->insert(pair<const Dof*, unsigned int>(dof[j], 0));
    }
  }
}

void DofManager::generateGlobalIdSpace(void){
  const map<const Dof*, unsigned int, DofComparator>::iterator
    end = globalId->end();

  map<const Dof*, unsigned int, DofComparator>::iterator
    it = globalId->begin();

  unsigned int id = 0;

  for(; it != end; it++){
    // Check if unknown
    if(it->second != isFixed){
      it->second = id;
      id++;
    }
  }

  serialize();
}

void DofManager::serialize(void){
  /*
  const map<const Dof*, unsigned int, DofComparator>::iterator end =
    --globalId->end();

  map<const Dof*, unsigned int, DofComparator>::iterator it =
    globalId->begin();

  unsigned int start = it->first->getEntity();
  unsigned int stop  = end->first->getEntity();

  cout << "Start: " << start            << endl
       << "Stop : " << stop             << endl
       << "Diff : " << stop - start + 1 << endl
       << "Tot  : " << globalId->size() << endl
       << endl;
  */
}

unsigned int DofManager::getGlobalId(const Dof& dof) const{
  const map<const Dof*, unsigned int, DofComparator>::
    iterator it = globalId->find(&dof);

  if(it == globalId->end())
    throw
      Exception("Their is no Dof %s", dof.toString().c_str());

  return it->second;
}

bool DofManager::fixValue(const Dof& dof, double value){
  // Get *REAL* Dof
  const map<const Dof*, unsigned int, DofComparator>::iterator it =
    globalId->find(&dof);

  // Check if 'dof' exists
  if(it == globalId->end())
    return false; // 'dof' doesn't exist

  // If 'dof' exists: it becomes fixed
  fixedDof->insert(pair<const Dof*, double>(it->first, value));
  it->second = isFixed;
  return true;
}

double DofManager::getValue(const Dof& dof) const{
  const map<const Dof*, unsigned int, DofComparator>::iterator end =
    globalId->end();

  map<const Dof*, unsigned int, DofComparator>::iterator it =
    globalId->find(&dof);

  if(it == end || it->second != isFixed)
    throw Exception("Unknown Dof: %s", dof.toString().c_str());

  else
    return fixedDof->find(it->first)->second;
}

string DofManager::toString(void) const{
  stringstream s;

  const map<const Dof*, unsigned int, DofComparator>::iterator end =
    globalId->end();

  map<const Dof*, unsigned int, DofComparator>::iterator it =
    globalId->begin();

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
