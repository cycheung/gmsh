#include "Mesh.h"
#include "Exception.h"

#include <list>
#include <sstream>

using namespace std;

Mesh::Mesh(const std::string fileName){ 
  // Alloc Memory //
  model       = new GModel("SmallFEM");
  physToGroup = new multimap<int, GroupOfElement*>;

  // Read Mesh //
  if(!model->readMSH(fileName))
    throw Exception("Can't open file: %s", fileName.c_str());

  // Get Entities (info) //
  vector<GEntity*> entity;
  model->getEntities(entity);
  nEntity = entity.size();

  // Get Entities (alloc) //
  group = new vector<GroupOfElement*>(nEntity);
  
  // Get Entities (get data) //
  for(unsigned int i = 0; i < nEntity; i++){
    (*group)[i] = new GroupOfElement(*(entity[i]), i);

    vector<int> physical = entity[i]->getPhysicalEntities();
    int nPhysical        = physical.size();

    for(int j = 0; j < nPhysical; j++)
      physToGroup->insert(pair<int, GroupOfElement*>(physical[j], (*group)[i]));
  }
}

Mesh::~Mesh(void){
  delete model;  
  delete physToGroup;

  for(unsigned int i = 0; i < nEntity; i++)
    delete (*group)[i];
  delete group;
}

const vector<GroupOfElement*> Mesh::getFromPhysical(int physical) const{
  pair<multimap<int, GroupOfElement*>::iterator, 
       multimap<int, GroupOfElement*>::iterator> startStop = 
    physToGroup->equal_range(physical);
  
  multimap<int, GroupOfElement*>::iterator it;
  list<GroupOfElement*> lst;
  
  for(it = startStop.first; it != startStop.second; it++)
    lst.push_back((*it).second);

  return vector<GroupOfElement*>(lst.begin(), lst.end());
}

string Mesh::toString(void) const{
  stringstream stream;

  stream << "*********************************************"    
	 << endl
	 << "*                    Mesh                   *"    
	 << endl
	 << "*********************************************"
	 << endl << endl;
  

  for(unsigned int i = 0; i < nEntity; i++)
    stream << (*group)[i]->toString() << endl;

  return stream.str();
}
