#include "Mesh.h"
#include "Exception.h"

#include <iostream>

using namespace std;

Mesh::Mesh(const std::string fileName){ 
  // Alloc Memory //
  model       = new GModel("SmallFEM");
  physToGroup = new map<int, Group*>;

  // Read Mesh //
  if(!model->readMSH(fileName))
    throw Exception("Can't open file: %s", fileName.c_str());

  // Get Entities //
  vector<GEntity*> entity;
  model->getEntities(entity);

  for(int i = 0; i < entity.size(); i++)
    cout << entity[i]->dim() << endl;
}
 
Mesh::~Mesh(void){
  delete model;
  delete physToGroup;
}

const vector<pair<int, Group*> > Mesh::getAllGroups(void) const{
  map<int, Group*>::iterator start = physToGroup->begin();
  map<int, Group*>::iterator stop  = physToGroup->end();  

  return std::vector<std::pair<int, Group*> >(start, stop);
}

string Mesh::toString(void) const{
  return string("Mesh");
}
