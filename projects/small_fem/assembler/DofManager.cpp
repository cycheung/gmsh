#include <sstream>
#include "DofManager.h"

using namespace std;
/*
DofManager::DofManager(const std::vector<MElement*>& element){
  // Init Lookup struct and GroupOfDof //
  nGroup = element.size();
  dof               = new vector<Dof*>(getNbDofFromElements(element));
  globalId          = new map<Dof*, int    , DofComparator>;
  dofToEntityLookup = new map<Dof*, Entity*, DofComparator>;
  group             = new vector<GroupOfDof*>(nGroup);
  physical          = new multimap<int, Dof*>;

  // Add Elements to DofManager //
  nextId = 0;

  dofLookup = new set<Dof*, DofComparator>;

  for(int i = 0; i < nGroup; i++)
    add(*element[i], i);

  delete dofLookup;

  nDof = dof->size();
}
*/

DofManager::DofManager(void){
}

DofManager::~DofManager(void){
  for(int i = 0; i < nGroup; i++)
    delete (*group)[i];
  delete group;

  delete globalId;
  //delete dofToEntityLookup;
  delete physical;

  for(int i = 0; i < nDof; i++)
    delete (*dof)[i];
  delete dof;
}
/*
void DofManager::add(MElement& element, int groupId){  
  const int type = element.getType();
  const int nEntity = element.nEntity();
  const std::vector<Entity*>& entity = element.getAllEntities();

  (*group)[groupId] = new GroupOfDof(nEntity, element.getId());

  for(int i = 0; i < nEntity; i++){
    pair<set<Dof*, DofComparator>::iterator, bool> p;
    Dof* tmp = new Dof(entity[i]->getId(), type);

    p = dofLookup->insert(tmp);
 
    if(p.second){
      (*dof)[nextId] = tmp;
      globalId->insert(pair<Dof*, int>(tmp, nextId));
      dofToEntityLookup->insert(pair<Dof*, Entity*>(tmp, entity[i]));
      
      (*group)[groupId]->add(tmp);

      if(entity[i]->gotPhysical())
	physical->insert(pair<int, Dof*>(entity[i]->getPhysical(), tmp));

      nextId += 1;
    }

    else{
      delete tmp; // Dof already exists
      (*group)[groupId]->add(*(p.first)); // Add real Dof
    }
  }
  
  (*group)[groupId]->jacobian(element);
  (*group)[groupId]->orientation(element.getAllOrientations());
}

int DofManager::getNbDofFromElements(const vector<MElement*>& element) const{
  set<int> entityLookup;
  const int N = element.size();
  
  for(int i = 0; i < N; i++){
    const vector<int>& id = element[i]->getAllEntitiesId();
    const int M = id.size();
    
    for(int j = 0; j < M; j++)
      entityLookup.insert(id[j]);
  }
    
  return entityLookup.size();
}
*/
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
