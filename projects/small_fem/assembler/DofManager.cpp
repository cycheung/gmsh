#include <sstream>
#include "Exception.h"
#include "GroupOfDof.h"

#include "DofManager.h"

using namespace std;

DofManager::DofManager(void){
  // Alloc Struct //
  globalId = new map<const Dof*, int, DofComparator>;
  fixedDof = new map<const Dof*, double, DofComparator>;

  // Dofs Numbering //
  nextId = 0;
}

DofManager::~DofManager(void){
  delete globalId;
  delete fixedDof;
}

void DofManager::addDof(const FunctionSpace& fs){
  // Get Groups //
  const vector<GroupOfDof*>& god = fs.getAllGroups();

  // Number Dof //
  const unsigned int nGoD = god.size();

  for(unsigned int i = 0; i < nGoD; i++){
    // Dof from god[i]
    const vector<const Dof*>& dof = god[i]->getAll();

    // Number Them
    const unsigned int nDof = dof.size();

    for(unsigned int j = 0; j < nDof; j++){
      pair<map<const Dof*, int, DofComparator>::iterator,
	   bool> newDof = 
	
	globalId->insert(pair<const Dof*, int>(dof[j], nextId));

      if(newDof.second)
	nextId++;
    }
    // NB: If the jth Dof was already in 'globalId', nothing append !
  }

}

int DofManager::getGlobalId(const Dof& dof) const{
  const map<const Dof*, int, DofComparator>::iterator it = 
    globalId->find(&dof);

  if(it == globalId->end())
    throw 
      Exception("Their is no Dof %s", dof.toString().c_str());

  else
    return it->second; 
}

bool DofManager::fixValue(const Dof& dof, double value){
  // Get *REAL* Dof
  map<const Dof*, int, DofComparator>::const_iterator it = 
    globalId->find(const_cast<Dof*>(&dof));

  // Check if 'dof' exists
  if(it == globalId->end())
    return false; // 'dof' doesn't exist

  // Insert *REAL* Dof 
  return fixedDof->insert(std::pair<const Dof*, double>(it->first, value)).second;
}

pair<bool, double> DofManager::getValue(const Dof& dof) const{
  map<const Dof*, double, DofComparator>::iterator end = 
    fixedDof->end();

  map<const Dof*, double, DofComparator>::iterator it = 
    fixedDof->find(&dof);

  if(it == end)
    return pair<bool, double>(false, 42);

  else
    return pair<bool, double>(true, it->second);
}

string DofManager::toString(void) const{
  stringstream s;
  map<const Dof*, int, DofComparator>::iterator i = 
    globalId->begin();
  
  map<const Dof*, int, DofComparator>::iterator end = 
    globalId->end();
  
  for(; i != end; i++)
    s << "("  << (*i).first->toString() 

      << ": " << (*i).second 
      << ")"  << endl;

  return s.str();
}
