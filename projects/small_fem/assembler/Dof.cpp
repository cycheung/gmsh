#include <sstream>

#include "Dof.h"
#include "Exception.h"

Dof::Dof(const unsigned int entity, const unsigned int type){
  this->entity  = entity;
  this->type    = type;
}

Dof::~Dof(void){
}

std::string Dof::toString(void) const{
  std::stringstream s;
  
  s << "(" << entity << ", " << type << ")";

  return s.str();
}
