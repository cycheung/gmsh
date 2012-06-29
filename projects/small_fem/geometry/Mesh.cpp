#include "Mesh.h"
#include "Exception.h"

using namespace std;

Mesh::Mesh(const std::string fileName){ 
  model = new GModel("SmallFEM");
  
  if(!model->readMSH(fileName))
    throw Exception("Can't open file: %s", fileName.c_str());

}
 
Mesh::~Mesh(void){
  delete model;
}

const vector<pair<int, Group*> > Mesh::getAllGroups(void) const{
  map<int, Group*>::iterator start = physToGroup->begin();
  map<int, Group*>::iterator stop  = physToGroup->end();  

  return std::vector<std::pair<int, Group*> >(start, stop);
}

string Mesh::toString(void) const{
  return string("Mesh");
}
