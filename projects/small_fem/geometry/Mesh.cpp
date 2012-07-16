#include "Mesh.h"
#include "GroupOfElements.h"
#include "Exception.h"

#include <sstream>

using namespace std;

Mesh::Mesh(const std::string fileName){ 
  // Alloc Memory //
  model       = new GModel("SmallFEM");
  physToGroup = new multimap<int, Group*>;

  // Read Mesh //
  if(!model->readMSH(fileName))
    throw Exception("Can't open file: %s", fileName.c_str());

  // Get Entities (info) //
  vector<GEntity*> entity;
  model->getEntities(entity);
  nEntity = entity.size();

  // Get Entities (alloc) //
  group = new vector<Group*>(nEntity);
  
  // Get Entities (get data) //
  for(unsigned int i = 0; i < nEntity; i++){
    (*group)[i] = new GroupOfElements(*(entity[i]), i);

    vector<int> physical = entity[i]->getPhysicalEntities();
    int nPhysical        = physical.size();

    for(int j = 0; j < nPhysical; j++)
      physToGroup->insert(pair<int, Group*>(physical[j], (*group)[i]));
  }
}
 
Mesh::~Mesh(void){
  delete model;  
  delete physToGroup;

  for(unsigned int i = 0; i < nEntity; i++)
    delete (*group)[i];
  delete group;
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
