#include <sstream>
#include "GroupOfElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"
#include "Exception.h"
#include "DofManager.h"

using namespace std;

DofManager::DofManager(FunctionSpace& fs){
  // Remember FunctionSpace & Associate//
  fs.associate(*this);
  this->fs = &fs;

  // Get Support from FunctionSpace //
  const GroupOfElement& support    = fs.getSupport();
  int nElement                     = support.getNumber();
  const vector<MElement*>& element = support.getAll();
  
  nTotVertex = support.getNVertex();

  // Init Struct //
  dof      = new set<Dof*, DofComparator>;         
  group    = new vector<GroupOfDof*>(nElement);
  globalId = new map<const Dof*, int, DofComparator>;
  fixedDof = new map<const Dof*, double, DofComparator>;

  elementToGroup = new map<const MElement*, 
			   const GroupOfDof*, 
			   ElementComparator>;

  // Create Dofs & Numbering//
  nextId = 0;

  // Loop over Element
  for(int i = 0; i < nElement; i++){
    // Get Dof for this Element
    vector<Dof*> myDof = fs.getKeys(*(element[i]));
    int nDof           = myDof.size();
    
    // Create new GroupOfDof
    GroupOfDof* god = new GroupOfDof(nDof, *(element[i])); 
    
    (*group)[i] = god;
    elementToGroup->insert
      (pair<const MElement*, const GroupOfDof*>(element[i], god));

    // Add Dof
    for(int j = 0; j < nDof; j++)
      insertDof(myDof[j], god);
  }
}

DofManager::~DofManager(void){
  int nElement = group->size();

  for(int i = 0; i < nElement; i++)
    delete (*group)[i];
  delete group;

  set<Dof*>::iterator it;
  set<Dof*>::iterator end = dof->end();

  for(it = dof->begin(); it != end; it++)
    delete *it;
  delete dof;

  delete globalId;
  delete fixedDof;
  delete elementToGroup;
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

const GroupOfDof& DofManager::getGroup(const MElement& element) const{
  const map<const MElement*, 
	    const GroupOfDof*, 
	    ElementComparator>::iterator it = 
    elementToGroup->find(&element);

  if(it == elementToGroup->end())
    throw 
      Exception("Their is no Element %d", element.getNum());

  else
    return *(it->second);   
}

void DofManager::insertDof(Dof* d, GroupOfDof* god){
  // Try to insert Dof //
  pair<set<Dof*, DofComparator>::iterator, bool> p;
  p = dof->insert(d);
 
  // If insertion is OK (Dof 'd' didn't exist) //
  //   --> Add new Dof
  if(p.second){
    globalId->insert(pair<Dof*, int>(d, nextId));
    
    god->add(d);
    
    nextId += 1;
  }
  
  // If insertion failed (Dof 'd' already exists) //
  //   --> delete 'd' and add existing Dof
  else{
    delete d; 
    god->add(*(p.first));
  }
}

bool DofManager::fixValue(const Dof& dof, double value){
  // Get *REAL* Dof
  set<Dof*, DofComparator>::const_iterator it = 
    this->dof->find(const_cast<Dof*>(&dof));

  // Check if 'dof' exists
  if(it == this->dof->end())
    return false; // 'dof' doesn't exist

  // Insert *REAL* Dof 
  return fixedDof->insert(std::pair<const Dof*, double>(*it, value)).second;
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
