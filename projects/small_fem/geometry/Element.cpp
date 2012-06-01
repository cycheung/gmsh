#include <sstream>
#include "Element.h"
#include "Exception.h"

Element::Element(const int id, const int type){
  this->id   = id;
  this->type = type;
  
  switch(type){
  case 30: numberOfEntity = 3; break;
  case 31: numberOfEntity = 3; break;
  default: throw Exception("Unknow ElementType"); break;
  }

  node        = new std::vector<Node*>(numberOfEntity);
  entity      = new std::vector<Entity*>(numberOfEntity);
  orientation = new std::vector<int>(numberOfEntity);
  entityId    = NULL;
}

Element::~Element(void){
  delete node;
  // Elements are not responsible for deleting nodes...
  delete entity;
  // Elements are not responsible for deleting entities...
  delete orientation;

  if(entityId)
    delete entityId;

  if(jac)
    delete jac;
}

const std::vector<int>& Element::getAllEntitiesId(void){
  if(!entityId){
    entityId = new std::vector<int>(numberOfEntity);

    for(int i = 0; i < numberOfEntity; i++)
      (*entityId)[i] = (*entity)[i]->getId();
  }
  
  return *entityId;
}

std::string Element::toString(void) const{
  std::stringstream stream;
  
  stream << "Element " << id << std::endl;
  stream << "******* " << std::endl;

  stream << "  -- Entity:" << std::endl;
  for(int i = 0; i < numberOfEntity; i++)
    stream << "    ++ " << (*entity)[i]->getId() 
	   << " ("      << (*entity)[i]->getType() << ")" 
	   << " ["      << (*orientation)[i] << "]"
	   << std::endl;

  return stream.str();
}

void Element::buildJacobian(void){
  jac = new Jacobian(*node);
}
