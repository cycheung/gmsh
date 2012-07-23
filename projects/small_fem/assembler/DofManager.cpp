#include <sstream>
#include "MVertex.h"
#include "DofManager.h"

using namespace std;

DofManager::DofManager(const GroupOfElement& goe){
  
  // Init Lookup struct and GroupOfDof //
  nGroup            = goe.getNumber();
  globalId          = new map<Dof*, int    , DofComparator>;
  group             = new vector<GroupOfDof*>(nGroup);

  // Get MElements //
  const vector<MElement*>& element = goe.getAll();

  // Add Elements to DofManager //
  nextId = 0;

  dofLookup = new set<Dof*, DofComparator>;

  for(int i = 0; i < nGroup; i++)
    add(*(element[i]), i);


  dof  = new vector<Dof*>(dofLookup->begin(), dofLookup->end());  
  nDof = dof->size();
  
  delete dofLookup;
}


DofManager::DofManager(void){
}

DofManager::~DofManager(void){
  for(int i = 0; i < nGroup; i++)
    delete (*group)[i];
  delete group;

  delete globalId;
  
  for(int i = 0; i < nDof; i++)
    delete (*dof)[i];
  delete dof;
}

void DofManager::add(MElement& element, int groupId){  
  // Up to now, we do only vertices ...//
  const int type = 1; 
  const int nEntity = element.getNumVertices();
  
  vector<MVertex*> entity;
  element.getVertices(entity);

  (*group)[groupId] = new GroupOfDof(nEntity, element);
  
  for(int i = 0; i < nEntity; i++){
    pair<set<Dof*, DofComparator>::iterator, bool> p;
    Dof* tmp = new Dof(entity[i]->getNum(), type);

    p = dofLookup->insert(tmp);
 
    if(p.second){
      globalId->insert(pair<Dof*, int>(tmp, nextId));
      
      (*group)[groupId]->add(tmp);

      nextId += 1;
    }

    else{
      delete tmp; // Dof already exists
      (*group)[groupId]->add(*(p.first)); // Add real Dof
    }
  }
  
}

string DofManager::toString(void) const{
  stringstream s;
  map<Dof*, int, DofComparator>::iterator i   = globalId->begin();
  map<Dof*, int, DofComparator>::iterator end = globalId->end();
  
  for(; i != end; i++)
    s << "("  << (*i).first->toString() 

      << ": " << (*i).second 
      << ")"  << endl;

  return s.str();
}
