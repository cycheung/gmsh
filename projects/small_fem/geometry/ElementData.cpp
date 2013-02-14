#include "ElementData.h"
#include "Exception.h"

ElementData::ElementData(void){
  hasGod = false;

  god = NULL;
}

ElementData::~ElementData(void){
}

void ElementData::setGroupOfDof(const GroupOfDof& god){
  this->hasGod = true;
  this->god    = &god;
}

const GroupOfDof& ElementData::getGroupOfDof(void) const{
  if(hasGod)
    return *god;

  else
    throw Exception("This ElementData isn't linked with a GroupOfDof");
}
